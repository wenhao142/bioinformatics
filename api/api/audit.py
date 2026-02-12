from datetime import datetime, timezone
from typing import List, Literal, TypedDict

from fastapi import APIRouter, Depends

from api.auth import current_user
from api.audit_store import get_store

EventType = Literal["upload", "run", "export"]


class AuditEvent(TypedDict):
    ts: str
    user: str
    project_id: str
    action: EventType
    detail: str


_EVENTS: List[AuditEvent] = []
_STORE = get_store()

router = APIRouter(prefix="/audit", tags=["audit"])


def log_event(user_email: str, project_id: str, action: EventType, detail: str = ""):
    evt = {
        "ts": datetime.now(timezone.utc).isoformat(),
        "user": user_email,
        "project_id": project_id,
        "action": action,
        "detail": detail,
    }
    _EVENTS.append(evt)
    if _STORE:
        _STORE.add_event(user_email, project_id, action, detail)


@router.get("")
def list_events(user=Depends(current_user)):
    if _STORE:
        events = _STORE.list_events(user["email"], user["role"] == "admin")
    else:
        events = _EVENTS if user["role"] == "admin" else [e for e in _EVENTS if e["user"] == user["email"]]
    return {"events": events}
