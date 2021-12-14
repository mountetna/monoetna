from flask_appbuilder.security.sqla.models import User


def test_fresh_login(admin_client):
  response = admin_client.get('/login/')
  assert response.status == '302 FOUND'
  assert response.headers['Location'] == 'https://janus.ucsf.edu/login?referer=http%3A%2F%2Flocalhost%2Flogin%2F'

# Simply validate the integration of all parts.
def test_with_cookie(admin_client, app, session):
  admin_client.set_cookie(
    'localhost',
    app.config['ETNA_AUTH_COOKIE_NAME'],
    'eyJhbGciOiJSUzI1NiJ9.eyJlbWFpbCI6InphY2hhcnkuY29sbGluc0B1Y3NmLmVkdSIsIm5hbWUiOiJaYWNoYXJ5IENvbGxpbnMiLCJwZXJtIjoiQTptdmlyMTthOmFkbWluaXN0cmF0aW9uLGNsdWVzMSxjb3Byb2plY3RzX3RlbXBsYXRlLHNpY2NhMSx0cmlhZ2UseGF1dDIseGNyczEseGhsdDEseGhsdDIseG5lbzEseG5ldTEseG9yZzEseHZpcjEsemFjaF90ZXN0X3Byb2plY3Rfb25lO2U6ZHNjb2xhYix0ZXN0X3Byb2plY3QiLCJmbGFncyI6ImluZ2VzdDt2dWxjYW47dGltdXJhZHZhbmNlZCIsImV4cCI6MTYzOTU1NzU1MX0.16BqHM2BCNDN0dPzFA6xMXe0H7tagaLlC7oPPjg-g6fHp3ehtk_dZh-BAsmNnqNXmo7_PyercjKQrCtbeAa-SHIdoLzqYu6nzvxrmRDk1p7tIhcncmVwAyjE5HBD740dPDiZzYjgnSx4M-cXULT1MOWzPZC2Eh87MvgtFRI7PCFYATm0nvH5otG9cm8iBtBKQkWpu5WIoYcP2jlLZ7bBSITWNKMC-jARA6iVhbquf02_129w9_7sab0uZ8NjvqjUIIWqyF1rAIul9c0ihtMrwvLEZB5LgknMqTe6IyVZfyepl0b5CUjy1f0d3vd9GQ6scm6VXxXTWjmnt_zgqD6Ze9UJFwsyT35gk7ECipbgw9fiBshMcHWE7sPcmeyh-fdahxz7JXiTB8PVkV0KUehY_oEPWVjUV2n0qhFeHvUsH3nziTQ5h2gVkXWL8ndQ6lr1Bt4_wDnV32OMsQNHg3ingAdutVALgE0QBZceoK_UaJ-JazDLCVLjKb8wzwU6gGO-NCL8uT21z2TD-yZBBgaLaYjvjQ7BKmbPH93ICacU2T6fxDYni3essDyoV4l2eo8wW4xTCyllrqXVWRGWV7efs3mbPkHHWLFDgvruXFu-Z4xAvm43wj8brFNEnN4yY6L0h1HR5mN-5KAKFcxIGDyIMA8ULypLhpOfdgHrWiDBDWc')
  response = admin_client.get('/login/')
  assert response.status == '302 FOUND'
  assert response.headers['Location'] == 'http://localhost/'

  users = session.query(User).order_by(User.id.desc()).limit(1).all()
  assert len(users) == 1
  user = users[0]

  assert 'Admin' in {r.name for r in user.roles}
  assert 'mvir1_viewer' in {r.name for r in user.roles}

