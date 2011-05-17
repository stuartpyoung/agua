CREATE TABLE project
(
    username    VARCHAR(20) NOT NULL,
    name        VARCHAR(20) NOT NULL,
    description VARCHAR(255) NOT NULL DEFAULT '',
    notes       TEXT NOT NULL DEFAULT '',

    PRIMARY KEY (username, name)
);
