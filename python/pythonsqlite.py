{\rtf1\ansi\ansicpg1252\cocoartf2511
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\froman\fcharset0 Times-Roman;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\sl280\partightenfactor0

\f0\fs24 \cf2 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 -- projects table\
CREATE TABLE IF NOT EXISTS projects (\
    id integer PRIMARY KEY,\
    name text NOT NULL,\
    begin_date text,\
    end_date text\
);\
 \
-- tasks table\
CREATE TABLE IF NOT EXISTS tasks (\
    id integer PRIMARY KEY,\
    name text NOT NULL,\
    priority integer,\
    project_id integer NOT NULL,\
    status_id integer NOT NULL,\
    begin_date text NOT NULL,\
    end_date text NOT NULL,\
    FOREIGN KEY (project_id) REFERENCES projects (id)\
);\
\uc0\u8232 \
}