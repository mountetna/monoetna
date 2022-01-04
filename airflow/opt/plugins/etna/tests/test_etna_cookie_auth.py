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
        "eyJhbGciOiJSUzI1NiJ9.eyJlbWFpbCI6ImRldmVsb3BlckB1Y3NmLmVkdSIsIm5hbWUiOiJEZXZlbG9wZXJQZXJzb24iLCJwZXJtIjoiYTphZG1pbmlzdHJhdGlvbixpcGksdGVzdC1wcm9qZWN0IiwiZmxhZ3MiOiJ2dWxjYW4iLCJleHAiOjI1MDM3MTM1NjZ9.1zO097iI_bNt4_Jn2pDP9VZxCGWGiBRFEz3TVjStY8spA_8tNjLWlhAzoJvWegpetnSmPVRDOCvNb6n8kjLv2LV716sdnoS1zuhtyMbbzLUHOpDPKpBuFVQN3fAiWjkw-Qo-4kQsMlk64pGK5GT4rQkUAouVlhsD1I7o-eVX9MKiTfwoHg5HndmjDFyRFsj6uTACAAYMGsgV3nGg3m6v4ebNeBM9AOYsBRnoR30Qd6E1j3VqxmD_RqPn57FQ-r6Ult-SfrK6lkldGxZa9xgAoxF-4s67C7AozcMvYq-A5ri4DbgwRgqwXf6kUNWpotGiUFPWqdKbRBLMH6e4TSZWuw",
    )
    response = admin_client.get("/login/")
    assert response.status == "302 FOUND"
    assert response.headers["Location"] == "http://localhost/"

    users = session.query(User).order_by(User.id.desc()).limit(1).all()
    assert len(users) == 1
    user = users[0]

    assert "Admin" in {r.name for r in user.roles}
    assert "ipi_viewer" in {r.name for r in user.roles}
