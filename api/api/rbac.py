from typing import Dict, Literal

from fastapi import APIRouter, Depends, HTTPException, Path, status
from pydantic import BaseModel

from api.auth import current_user
from api.audit import log_event

Role = Literal["admin", "analyst", "viewer"]

# In-memory store for demo; replace with DB later.
PROJECTS: Dict[str, Dict[str, Role]] = {}

router = APIRouter(prefix="/projects", tags=["projects"])


class ProjectCreate(BaseModel):
    project_id: str


class MemberAdd(BaseModel):
    email: str
    role: Role


def make_role_dependency(allowed: tuple[Role, ...]):
    def dep(project_id: str = Path(..., description="Project identifier"), user=Depends(current_user)):
        members = PROJECTS.get(project_id, {})
        role = members.get(user["email"])
        if role not in allowed:
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Forbidden for your role")
        return user

    return dep


@router.post("", status_code=201)
def create_project(body: ProjectCreate, user=Depends(current_user)):
    if user["role"] != "admin":
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Admin only")
    PROJECTS.setdefault(body.project_id, {})[user["email"]] = "admin"
    return {"project_id": body.project_id}


@router.post("/{project_id}/members", status_code=201)
def add_member(
    body: MemberAdd,
    project_id: str = Path(..., description="Project identifier"),
    user=Depends(current_user),
):
    if user["role"] != "admin":
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Admin only")
    if project_id not in PROJECTS:
        raise HTTPException(status_code=404, detail="Project not found")
    PROJECTS[project_id][body.email] = body.role
    return {"project_id": project_id, "email": body.email, "role": body.role}


@router.get("/{project_id}/me")
def my_role(project_id: str, user=Depends(current_user)):
    role = PROJECTS.get(project_id, {}).get(user["email"])
    if role is None:
        raise HTTPException(status_code=404, detail="Not a member")
    return {"project_id": project_id, "email": user["email"], "role": role}


upload_dep = make_role_dependency(("admin", "analyst"))
view_dep = make_role_dependency(("admin", "analyst", "viewer"))


@router.post("/{project_id}/upload")
def upload_stub(project_id: str, user=Depends(upload_dep)):
    log_event(user["email"], project_id, "upload", detail="upload stub")
    return {"project_id": project_id, "action": "upload", "by": user["email"]}


@router.post("/{project_id}/run")
def run_stub(project_id: str, user=Depends(upload_dep)):
    log_event(user["email"], project_id, "run", detail="run stub")
    return {"project_id": project_id, "action": "run", "by": user["email"]}


@router.get("/{project_id}/view")
def view_stub(project_id: str, user=Depends(view_dep)):
    return {"project_id": project_id, "action": "view", "by": user["email"]}


@router.post("/{project_id}/export")
def export_stub(project_id: str, user=Depends(view_dep)):
    log_event(user["email"], project_id, "export", detail="export stub")
    return {"project_id": project_id, "action": "export", "by": user["email"]}
