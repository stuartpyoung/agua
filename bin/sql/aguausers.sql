CREATE TABLE IF NOT EXISTS aguausers
(
    username        VARCHAR(30) NOT NULL,
    password        VARCHAR(30) NOT NULL,
    groupname       VARCHAR(50),
    email           VARCHAR(50),
    firstname       VARCHAR(20),
    lastname        VARCHAR(20),
    description     TEXT,
    homedir         VARCHAR(255) NOT NULL,
    datetime        DATETIME NOT NULL,

    PRIMARY KEY (username)
);