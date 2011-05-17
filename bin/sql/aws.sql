CREATE TABLE IF NOT EXISTS aws
(
    username        VARCHAR(30) NOT NULL,
    volumesize      INT(6),
    amazonuserid    VARCHAR(30) NOT NULL,
    ec2publiccert   TEXT,
    ec2privatekey   TEXT,
    awsaccesskeyid  VARCHAR(100),
    secretaccesskey VARCHAR(100),

    PRIMARY KEY (username)
);