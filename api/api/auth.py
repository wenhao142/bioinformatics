import os
import time
import hashlib
import hmac
import secrets
from typing import Any

import jwt
import psycopg2
from fastapi import APIRouter, Depends, HTTPException, Security, status
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from pydantic import BaseModel

router = APIRouter(prefix="/auth", tags=["auth"])
security = HTTPBearer()
TRUE_VALUES = {"1", "true", "yes", "on"}

JWT_SECRET = os.getenv("JWT_SECRET", "dev-secret-change-me")
JWT_EXPIRE_SECONDS = int(os.getenv("JWT_EXPIRE_SECONDS", "3600"))
ADMIN_EMAIL = os.getenv("ADMIN_EMAIL", "admin@example.com")
ADMIN_PASSWORD = os.getenv("ADMIN_PASSWORD", "password")
DEFAULT_USER_PASSWORD = os.getenv("DEFAULT_USER_PASSWORD", "password")
PASSWORD_ROUNDS = int(os.getenv("PASSWORD_HASH_ROUNDS", "240000"))
VALID_ROLES = {"admin", "analyst", "viewer"}


def _hash_password(password: str, salt: str | None = None) -> str:
    if not password:
        raise ValueError("Password cannot be empty")
    final_salt = salt or secrets.token_hex(16)
    digest = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), final_salt.encode("utf-8"), PASSWORD_ROUNDS)
    return f"pbkdf2_sha256${PASSWORD_ROUNDS}${final_salt}${digest.hex()}"


def _verify_password(password: str, stored_hash: str) -> bool:
    parts = stored_hash.split("$")
    if len(parts) != 4 or parts[0] != "pbkdf2_sha256":
        # Backward-compatible fallback for plain-text legacy values.
        return hmac.compare_digest(password, stored_hash)
    _, rounds_raw, salt, digest_hex = parts
    try:
        rounds = int(rounds_raw)
    except ValueError:
        return False
    candidate = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), salt.encode("utf-8"), rounds).hex()
    return hmac.compare_digest(candidate, digest_hex)


class MemoryUserStore:
    def __init__(self):
        self.users: dict[str, dict[str, str]] = {}

    def get_user(self, email: str) -> dict[str, str] | None:
        row = self.users.get(email.lower())
        return dict(row) if row else None

    def create_user(self, email: str, password_hash: str, role: str) -> bool:
        key = email.lower()
        if key in self.users:
            return False
        self.users[key] = {"email": email, "password_hash": password_hash, "role": role}
        return True

    def upsert_user(self, email: str, password_hash: str, role: str) -> None:
        self.users[email.lower()] = {"email": email, "password_hash": password_hash, "role": role}

    def list_users(self) -> list[dict[str, str]]:
        rows = [{"email": row["email"], "role": row["role"]} for row in self.users.values()]
        rows.sort(key=lambda item: item["email"])
        return rows


class PgUserStore:
    def __init__(self):
        self.conn = psycopg2.connect(self._dsn())
        self.conn.autocommit = True
        self._ensure_table()

    def _dsn(self) -> str:
        dsn = os.getenv("DATABASE_URL")
        if dsn:
            return dsn
        return (
            f"dbname={os.getenv('POSTGRES_DB', 'omics')} "
            f"user={os.getenv('POSTGRES_USER', 'postgres')} "
            f"password={os.getenv('POSTGRES_PASSWORD', 'postgres')} "
            f"host={os.getenv('POSTGRES_HOST', 'localhost')} "
            f"port={os.getenv('POSTGRES_PORT', '5432')}"
        )

    def _ensure_table(self) -> None:
        with self.conn.cursor() as cur:
            cur.execute(
                """
                CREATE TABLE IF NOT EXISTS app_users (
                    email TEXT PRIMARY KEY,
                    password_hash TEXT NOT NULL,
                    role TEXT NOT NULL,
                    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
                    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
                );
                CREATE INDEX IF NOT EXISTS app_users_role_idx ON app_users(role);
                """
            )

    def get_user(self, email: str) -> dict[str, str] | None:
        with self.conn.cursor() as cur:
            cur.execute("SELECT email, password_hash, role FROM app_users WHERE lower(email)=lower(%s)", (email,))
            row = cur.fetchone()
        if not row:
            return None
        return {"email": row[0], "password_hash": row[1], "role": row[2]}

    def create_user(self, email: str, password_hash: str, role: str) -> bool:
        with self.conn.cursor() as cur:
            cur.execute(
                """
                INSERT INTO app_users(email, password_hash, role)
                VALUES (%s, %s, %s)
                ON CONFLICT (email) DO NOTHING
                RETURNING email
                """,
                (email, password_hash, role),
            )
            created = cur.fetchone()
        return bool(created)

    def upsert_user(self, email: str, password_hash: str, role: str) -> None:
        with self.conn.cursor() as cur:
            cur.execute(
                """
                INSERT INTO app_users(email, password_hash, role)
                VALUES (%s, %s, %s)
                ON CONFLICT (email) DO UPDATE
                SET password_hash=EXCLUDED.password_hash,
                    role=EXCLUDED.role,
                    updated_at=NOW()
                """,
                (email, password_hash, role),
            )

    def list_users(self) -> list[dict[str, str]]:
        with self.conn.cursor() as cur:
            cur.execute("SELECT email, role FROM app_users ORDER BY email")
            rows = cur.fetchall()
        return [{"email": row[0], "role": row[1]} for row in rows]


def _should_use_db() -> bool:
    return os.getenv("AUTH_USE_DB", "false").lower() in TRUE_VALUES


def _seed_default_users(store: Any) -> None:
    defaults = [
        (ADMIN_EMAIL, ADMIN_PASSWORD, "admin"),
        ("analyst@example.com", DEFAULT_USER_PASSWORD, "analyst"),
        ("viewer@example.com", DEFAULT_USER_PASSWORD, "viewer"),
    ]
    for email, password, role in defaults:
        existing = store.get_user(email)
        if existing:
            continue
        store.create_user(email=email, password_hash=_hash_password(password), role=role)


def _create_store():
    if _should_use_db():
        try:
            store = PgUserStore()
            _seed_default_users(store)
            return store
        except Exception:
            pass
    store = MemoryUserStore()
    _seed_default_users(store)
    return store


USER_STORE = _create_store()


class LoginRequest(BaseModel):
    email: str
    password: str


class RegisterRequest(BaseModel):
    email: str
    password: str
    role: str = "viewer"


class TokenResponse(BaseModel):
    access_token: str
    token_type: str = "bearer"
    expires_in: int


def create_token(email: str, role: str) -> str:
    now = int(time.time())
    payload = {
        "sub": email,
        "role": role,
        "iat": now,
        "exp": now + JWT_EXPIRE_SECONDS,
    }
    return jwt.encode(payload, JWT_SECRET, algorithm="HS256")


@router.post("/login", response_model=TokenResponse)
def login(req: LoginRequest):
    user = USER_STORE.get_user(req.email)
    if not user or not _verify_password(req.password, user["password_hash"]):
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid credentials")
    token = create_token(user["email"], user["role"])
    return TokenResponse(access_token=token, expires_in=JWT_EXPIRE_SECONDS)


def current_user(credentials: HTTPAuthorizationCredentials = Security(security)) -> dict:
    token = credentials.credentials
    try:
        payload = jwt.decode(token, JWT_SECRET, algorithms=["HS256"])
    except jwt.PyJWTError:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid token")
    email = payload.get("sub")
    role = payload.get("role")
    if not isinstance(email, str) or not isinstance(role, str) or role not in VALID_ROLES:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid token payload")
    return {"email": email, "role": role}


@router.get("/me")
def me(user=Depends(current_user)):
    return {"email": user["email"], "role": user["role"]}


@router.post("/register", status_code=201)
def register(req: RegisterRequest, user=Depends(current_user)):
    if user["role"] != "admin":
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Admin only")

    email = req.email.strip()
    role = req.role.strip().lower()
    if not email:
        raise HTTPException(status_code=400, detail="Email is required")
    if role not in VALID_ROLES:
        raise HTTPException(status_code=400, detail=f"Role must be one of: {', '.join(sorted(VALID_ROLES))}")
    if len(req.password) < 8:
        raise HTTPException(status_code=400, detail="Password must be at least 8 characters")

    created = USER_STORE.create_user(email, _hash_password(req.password), role)
    if not created:
        raise HTTPException(status_code=409, detail="User already exists")
    return {"email": email, "role": role}


@router.get("/users")
def list_users(user=Depends(current_user)):
    if user["role"] != "admin":
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Admin only")
    return {"users": USER_STORE.list_users()}


def reset_auth_store_for_tests() -> None:
    global USER_STORE
    USER_STORE = MemoryUserStore()
    _seed_default_users(USER_STORE)
