CREATE TABLE IF NOT EXISTS view
(
    username        VARCHAR(30) NOT NULL,
    project         VARCHAR(20) NOT NULL,
    view            VARCHAR(20) NOT NULL,
    species         VARCHAR(20),
    build           VARCHAR(20),
    chromosome      VARCHAR(20),
    tracklist       TEXT,
    start           INT(15),
    stop            INT(15),
    notes           TEXT,
    datetime        DATETIME NOT NULL,

    PRIMARY KEY (username, project, view)
);
