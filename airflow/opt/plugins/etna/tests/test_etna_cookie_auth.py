from flask_appbuilder.security.sqla.models import User


def test_fresh_login(admin_client):
    response = admin_client.get("/login/")
    assert response.status == "302 FOUND"
    assert (
        response.headers["Location"]
        == "https://janus.development.local/login?referer=http%3A%2F%2Flocalhost%2Flogin%2F"
    )


# Simply validate the integration of all parts.
def test_with_cookie(admin_client, app, session):
    admin_client.set_cookie(
        "localhost",
        app.config["ETNA_AUTH_COOKIE_NAME"],
        "eyJhbGciOiJSUzI1NiJ9.eyJlbWFpbCI6ImRldmVsb3BlckB1Y3NmLmVkdSIsIm5hbWUiOiJEZXZlbG9wZXJQZXJzb24iLCJwZXJtIjoiYTphZG1pbmlzdHJhdGlvbixpcGksdGVzdC1wcm9qZWN0IiwiZmxhZ3MiOiJ2dWxjYW4iLCJleHAiOjE2Mzk2MTE2MTB9.1WH9ld8OocS2il-HHp2-FhtYYkgmHtRvGJ0pppHU7MKnkZnrmJbrrgfnNQXqnD7ng8iZnrOvDYPetK9nnMb3WU2lF8CAalw9osHyw5hnvZLPhdcZyuDq18Fe2SYwAEmKAAytN_n4N8D5yIjRGygBH-6V9whZ4bzx2jUrchIquUgugAeE2Ux9qn5YfM359GEtJ1xizx4GtzImxtfOAmqt12oWFb_FwunFWLzusXPumk2-t1beyZclH9-M3OYmtFNCy9fg79lyb8slqlgUw4LdqrVtq4hKDwV-2C-RttW2bMUuQnuDptoBeL8b_DD4iZ4Mr06UG7XvWlCkNL7sgQ5YsA",
    )
    response = admin_client.get("/login/")
    assert response.status == "302 FOUND"
    assert response.headers["Location"] == "http://localhost/"

    users = session.query(User).order_by(User.id.desc()).limit(1).all()
    assert len(users) == 1
    user = users[0]

    assert "Admin" in {r.name for r in user.roles}
    assert "ipi_viewer" in {r.name for r in user.roles}
