CREATE TABLE IF NOT EXISTS workflow
(
    username            VARCHAR(30) NOT NULL,
    project             VARCHAR(20) NOT NULL,
    name                VARCHAR(20) NOT NULL,
    number              INT(12),

    description         TEXT,
    notes               TEXT,

    PRIMARY KEY  (username, project, name, number)
);