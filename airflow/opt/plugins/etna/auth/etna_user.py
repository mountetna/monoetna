import dataclasses
import re
from typing import Dict, List, Callable, Optional

from etna.utils.parser import (
    parse_many_joined_by,
    Parser,
    literal_parser,
    take_left,
    eof_parser,
    regex_parser,
    parser_apply,
    parser_lift,
    either_parser,
)

ROLE_MAP = dict(
    a="admin",
    v="viewer",
    e="editor",
)


@dataclasses.dataclass
class Permission:
    project_name: str
    see_restricted: bool
    role: str


@dataclasses.dataclass
class EtnaUser:
    name: str
    email: str
    flags: List[str]
    perm: List[Permission]

    @property
    def username(self):
        return self.email.split("@")[0]

    @staticmethod
    def from_dict(payload: Dict[str, str]) -> "EtnaUser":
        name = payload["name"]
        email = payload["email"]
        flags = payload.get("flags", "").split(";")
        perm = parse_permissions(payload["perm"])

        return EtnaUser(
            name=name,
            email=email,
            flags=flags,
            perm=perm,
        )

    @property
    def projects(self) -> List[str]:
        return [p.project_name for p in self.perm]

    def project_permission(self, project_name: str) -> Optional[Permission]:
        for p in self.perm:
            if p.project_name == project_name:
                return p

    def has_any_role(self, project_name: str, *roles: str) -> bool:
        permission = self.project_permission(project_name)
        if permission is None:
            return False

        return permission.role in roles

    @property
    def is_superuser(self) -> bool:
        return self.has_any_role("administration", "admin")

    @property
    def is_supereditor(self) -> bool:
        return self.has_any_role("administration", "admin", "editor")

    @property
    def is_superviewer(self) -> bool:
        return self.has_any_role("administration", "admin", "editor", "viewer")

    def can_edit(self, project_name) -> bool:
        return self.is_supereditor or self.has_any_role(project_name, "admin", "editor")

    def can_view(self, project_name) -> bool:
        return self.is_superviewer or self.has_any_role(
            project_name, "admin", "editor", "viewer"
        )

    def is_admin(self, project_name) -> bool:
        return self.is_superuser or self.has_any_role(project_name, "admin")


def construct_permission(role: str) -> Callable[[List[str]], List[Permission]]:
    def construct_permissions(project_names: List[str]) -> List[Permission]:
        return [
            Permission(
                role=ROLE_MAP[role.lower()],
                see_restricted=role == role.upper(),
                project_name=project_name,
            )
            for project_name in project_names
        ]

    return construct_permissions


role_parser: Parser[str] = take_left(
    regex_parser(re.compile(r"[aAeEvV]")), literal_parser(":")
)
permission_parser: Parser[List[Permission]] = parser_apply(
    parser_apply(parser_lift(construct_permission), role_parser),
    parse_many_joined_by(regex_parser(re.compile(r"[^,\s;]+")), literal_parser(",")),
)


permissions_parser: Parser[List[List[Permission]]] = take_left(
    parse_many_joined_by(permission_parser, literal_parser(";")), eof_parser
)


def parse_permissions(perm: str) -> List[Permission]:
    result, error = permissions_parser.run_parser(perm)

    if result:
        return [p for ps in result for p in ps]

    raise error
