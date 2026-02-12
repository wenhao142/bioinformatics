import os
import time
from typing import Optional

import jwt
from fastapi import APIRouter, Depends, HTTPException, Security, status
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from pydantic import BaseModel

router = APIRouter(prefix="/auth", tags=["auth"])
security = HTTPBearer()

JWT_SECRET = os.getenv("JWT_SECRET", "dev-secret-change-me")
JWT_EXPIRE_SECONDS = int(os.getenv("JWT_EXPIRE_SECONDS", "3600"))
ADMIN_EMAIL = os.getenv("ADMIN_EMAIL", "admin@example.com")
ADMIN_PASSWORD = os.getenv("ADMIN_PASSWORD", "password")
DEFAULT_USER_PASSWORD = os.getenv("DEFAULT_USER_PASSWORD", "password")

# In-memory demo user store (replace with DB later)
DEMO_USERS = {
    ADMIN_EMAIL: {
        "password": ADMIN_PASSWORD,
        "role": "admin",
    }
}
for demo_email, role in (
    ("analyst@example.com", "analyst"),
    ("viewer@example.com", "viewer"),
):
    DEMO_USERS.setdefault(demo_email, {"password": DEFAULT_USER_PASSWORD, "role": role})


class LoginRequest(BaseModel):
    email: str
    password: str


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
    user = DEMO_USERS.get(req.email)
    if not user or user["password"] != req.password:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid credentials")
    token = create_token(req.email, user["role"])
    return TokenResponse(access_token=token, expires_in=JWT_EXPIRE_SECONDS)


def current_user(credentials: HTTPAuthorizationCredentials = Security(security)) -> dict:
    token = credentials.credentials
    try:
        payload = jwt.decode(token, JWT_SECRET, algorithms=["HS256"])
    except jwt.PyJWTError:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid token")
    return {"email": payload.get("sub"), "role": payload.get("role")}


@router.get("/me")
def me(user=Depends(current_user)):
    return {"email": user["email"], "role": user["role"]}
