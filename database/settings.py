from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker

DATABASE = 'sqlite:///db.sqlite3'
Engine = create_engine(
  DATABASE,
  echo=False
)
Base = declarative_base()
session = scoped_session(
    sessionmaker(
        autocommit = False,
    	autoflush = False,
    	bind = Engine
    )
)
Base = declarative_base()
Base.query = session.query_property()
