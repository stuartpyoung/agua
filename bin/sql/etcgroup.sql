CREATE TABLE IF NOT EXISTS etcgroup
(
    username        VARCHAR(30) NOT NULL,
    gid                 INT(9),
    groupname            VARCHAR(255),
    email               VARCHAR(50),
    firstname           VARCHAR(20),
    lastname            VARCHAR(20),
    description         TEXT,

    PRIMARY KEY (username, gid, groupname)
);
    