#!venv/bin/python
from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base, sessionmaker
from sqlalchemy_utils import database_exists, create_database

import settings

Base = declarative_base()

def connect_database():
    engine = create_engine(settings.database)
    if not database_exists(engine.url):
        print("## Making new database ##")
        create_database(engine.url)
        Base.metadata.create_all(engine)

    return sessionmaker(bind=engine)
