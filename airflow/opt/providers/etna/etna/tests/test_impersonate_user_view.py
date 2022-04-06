def test_form_shows_username_choices(app, admin_client):
    user = app.appbuilder.sm.find_user(username='test_admin')
    with admin_client.session_transaction() as sess:
        sess['user_id'] = user.id
        sess['_fresh'] = True  # https://flask-login.readthedocs.org/en/latest/#fresh-logins

    response = admin_client.get("/impersonateuserview/form")

    assert 'test_project_viewer_first_name' in response.data.decode('utf-8')


def test_form_submit_changes_user_roles(app, admin_client, reset_admin_role):
    with admin_client.session_transaction() as sess:
        admin_user = app.appbuilder.sm.find_user(username='test_admin')
        project_user = app.appbuilder.sm.find_user(username='test_project_viewer')

        admin_roles = admin_user.roles

        assert admin_user.roles != project_user.roles

        sess['user_id'] = admin_user.id
        sess['_fresh'] = True  # https://flask-login.readthedocs.org/en/latest/#fresh-logins

    response = admin_client.post("/impersonateuserview/form", data=dict(
        email=project_user.email
    ))
    assert response.status_code == 302 # redirect to home page after submit

    admin_user = app.appbuilder.sm.find_user(username='test_admin')
    project_user = app.appbuilder.sm.find_user(username='test_project_viewer')

    assert admin_user.roles == project_user.roles


def test_admin_user_can_access_view(app, admin_client):
    user = app.appbuilder.sm.find_user(username='test_admin')
    with admin_client.session_transaction() as sess:
        sess['user_id'] = user.id
        sess['_fresh'] = True  # https://flask-login.readthedocs.org/en/latest/#fresh-logins

    response = admin_client.get("/impersonateuserview/form")
    assert response.status_code == 200


def test_project_user_cannot_access_view(app, admin_client):
    with admin_client.session_transaction() as sess:
        sess['user_id'] = app.appbuilder.sm.find_user(username='test_project_viewer').id
        sess['_fresh'] = True  # https://flask-login.readthedocs.org/en/latest/#fresh-logins

    response = admin_client.get("/impersonateuserview/form")
    assert response.status_code == 302