from etna.auth.etna_user import (
    parse_permissions,
    permission_parser,
    EtnaUser,
    Permission,
)
from etna.utils.parser import take_left, eof_parser


def test_user_from_payload():
    payload = {
        "email": "zachary.collins@ucsf.edu",
        "name": "Zachary Collins",
        "perm": "A:mvir1",
        "flags": "ingest;vulcan;timuradvanced",
        "exp": 1639557551,
    }

    assert EtnaUser.from_dict(payload) == EtnaUser(
        name="Zachary Collins",
        email="zachary.collins@ucsf.edu",
        perm=[Permission(project_name="mvir1", see_restricted=True, role="admin")],
        flags=["ingest", "vulcan", "timuradvanced"],
    )


def test_parse_permissions():
    perm = "A:mvir1;a:administration,clues1,coprojects_template,sicca1,triage,xaut2,xcrs1,xhlt1,xhlt2,xneo1,xneu1,xorg1,xvir1,zach_test_project_one;e:dscolab,test_project"
    permissions = parse_permissions(perm)

    assert permissions == [
        Permission(project_name="mvir1", see_restricted=True, role="admin"),
        Permission(project_name="administration", see_restricted=False, role="admin"),
        Permission(project_name="clues1", see_restricted=False, role="admin"),
        Permission(
            project_name="coprojects_template", see_restricted=False, role="admin"
        ),
        Permission(project_name="sicca1", see_restricted=False, role="admin"),
        Permission(project_name="triage", see_restricted=False, role="admin"),
        Permission(project_name="xaut2", see_restricted=False, role="admin"),
        Permission(project_name="xcrs1", see_restricted=False, role="admin"),
        Permission(project_name="xhlt1", see_restricted=False, role="admin"),
        Permission(project_name="xhlt2", see_restricted=False, role="admin"),
        Permission(project_name="xneo1", see_restricted=False, role="admin"),
        Permission(project_name="xneu1", see_restricted=False, role="admin"),
        Permission(project_name="xorg1", see_restricted=False, role="admin"),
        Permission(project_name="xvir1", see_restricted=False, role="admin"),
        Permission(
            project_name="zach_test_project_one", see_restricted=False, role="admin"
        ),
        Permission(project_name="dscolab", see_restricted=False, role="editor"),
        Permission(project_name="test_project", see_restricted=False, role="editor"),
    ]


def test_parse_permission():
    permissions = take_left(permission_parser, eof_parser).run_parser(
        "A:project1,project2"
    )[0]
    assert permissions == [
        Permission(project_name="project1", see_restricted=True, role="admin"),
        Permission(project_name="project2", see_restricted=True, role="admin"),
    ]

    permissions = take_left(permission_parser, eof_parser).run_parser("e:project3")[0]
    assert permissions == [
        Permission(project_name="project3", see_restricted=False, role="editor"),
    ]

    permissions = take_left(permission_parser, eof_parser).run_parser(
        "v:project5,project1,project2"
    )[0]
    assert permissions == [
        Permission(project_name="project5", see_restricted=False, role="viewer"),
        Permission(project_name="project1", see_restricted=False, role="viewer"),
        Permission(project_name="project2", see_restricted=False, role="viewer"),
    ]
