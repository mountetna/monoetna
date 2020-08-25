-- BIG WARNING! MAKE SURE YOU ROTATE THE KEYS BELOW!
-- I am only including them here as an example.
-- These keys are not in production and never should be!

INSERT INTO private.apps (app_key, app_name)
VALUES
('BLAH_BLAH_qYk5HgHkq9ZXL7p8xwKXuYYNgwRt2tvuwsJMbjcm', 'janus'),
('yyQo35xKP4BLAH_BLAH_WRZxxlg1dFK6d74JxHMRDb8iWmDCwR', 'metis'),
('jsY1SKg9C3fbjL44lHF2BLAH_BLAH_PoeQO7g4LMWHypXb2Hnn', 'magma'),
('ydCVgJf6xoZ0wTZVIQKdR2ZJ0Kwu3vBLAH_BLAH_L7BG4s72JW', 'timur'),
('3OFChKR0rOUvMeAnUPYL9Wqzz2UOR81ONS6yr57WBLAH_BLAH_', 'polyphemus');

INSERT INTO private.groups (group_name)
VALUES
('administration'),
('some initial group'),
('another group');

INSERT INTO private.projects (group_id, project_name_full, project_name, project_description)
VALUES
(1, 'administration', 'Administration', 'The main project for Mount Etna.'),
(2, 'some special project', 'ssp', 'Some project of yours.'),
(3, 'another neat thing', 'ant', 'Another project of mine.');

INSERT INTO private.users (email, first_name, last_name)
VALUES
('some.guy@ucsf.edu', 'Hello', 'All'),
('some.gal@ucsf.edu', 'You', 'People');

UPDATE private.users
SET pass_hash = 'BLAH_BLAH_p8FLdIuPfZfqNYA1mpsuMQcYJXgTyVsG1Q89DT4A'
WHERE email = 'some.guy@ucsf.edu';

UPDATE private.users
SET pass_hash = 'BLAH_BLAH_mxq5HmykLwo6M3nR7rt9p8vBKYztsJHFgBvL0Arl'
WHERE email = 'some.gal@ucsf.edu';

INSERT INTO private.permissions (user_id, project_id, role)
VALUES
(1, 1, 'administrator'),
(2, 1, 'administrator'),
(1, 2, 'administrator'),
(2, 2, 'administrator'),
(1, 3, 'administrator'),
(2, 3, 'administrator');
