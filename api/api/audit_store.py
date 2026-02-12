import os
import psycopg2
from typing import Optional


class PgAuditStore:
    def __init__(self):
        self.conn = None
        self._connect()
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

    def _connect(self):
        self.conn = psycopg2.connect(self._dsn())
        self.conn.autocommit = True

    def _ensure_table(self):
        with self.conn.cursor() as cur:
            cur.execute(
                """
                CREATE TABLE IF NOT EXISTS audit_events (
                    id SERIAL PRIMARY KEY,
                    ts TIMESTAMPTZ NOT NULL DEFAULT now(),
                    user_email TEXT NOT NULL,
                    project_id TEXT NOT NULL,
                    action TEXT NOT NULL,
                    detail TEXT NOT NULL
                );
                """
            )

    def add_event(self, user_email: str, project_id: str, action: str, detail: str):
        with self.conn.cursor() as cur:
            cur.execute(
                "INSERT INTO audit_events (user_email, project_id, action, detail) VALUES (%s, %s, %s, %s)",
                (user_email, project_id, action, detail),
            )

    def list_events(self, user_email: Optional[str], is_admin: bool):
        with self.conn.cursor() as cur:
            if is_admin:
                cur.execute(
                    "SELECT ts, user_email, project_id, action, detail FROM audit_events ORDER BY ts DESC LIMIT 1000"
                )
            else:
                cur.execute(
                    "SELECT ts, user_email, project_id, action, detail FROM audit_events WHERE user_email=%s ORDER BY ts DESC LIMIT 1000",
                    (user_email,),
                )
            return [
                {
                    "ts": str(row[0]),
                    "user": row[1],
                    "project_id": row[2],
                    "action": row[3],
                    "detail": row[4],
                }
                for row in cur.fetchall()
            ]


def get_store():
    if os.getenv("AUDIT_USE_DB", "false").lower() in ("1", "true", "yes", "on"):
        try:
            return PgAuditStore()
        except Exception:
            # Fallback to None -> memory store
            return None
    return None
