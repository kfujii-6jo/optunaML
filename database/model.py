from datetime import datetime
from sqlalchemy import create_engine, Column
from sqlalchemy.types import Integer, Float, String, DateTime
# from sqlalchemy.dialects.mysql import INTEGER as Integer
# from sqlalchemy.dialects.mysql import TINYINT as Tinyint
# from sqlalchemy.dialects.mysql import TIMESTAMP as Timestamp
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker, relationship
from sqlalchemy.schema import ForeignKey

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

class Lattice(Base):
    __tablename__ = 'lattice'

    id = Column('id', Integer, primary_key=True, autoincrement=True, nullable=False)
    lattice_name = Column('lattice_name', String)
    cif_filepath = Column('cif_filepath', String)
    mol_filepath = Column('mol_filepath', String) 
    molecule_fingerprint = Column('molecule_fingerprint', String)
    A = Column('A', Float)
    B = Column('B', Float)
    C = Column('C', Float)
    dihedral_angle = Column('dihedral_angle', Float)

    def __repr__(self):
        return "<Lattice('id={}', lattice_name={}, cif_filepath={}, mol_filepath={}, molecule_fingerprint={}, A={}, B={}, C={}, dihedral_angle={} )>".format(
            self.id,
            self.lattice_name,
            self.cif_filepath,
            self.mol_filepath,
            self.molecule_fingerprint,
            self.A,
            self.B,
            self.C,
            self.dihedral_angle
        )

class FingerprintType(Base):
    __tablename__ = 'fingerprint_type'
    
    id = Column('id', Integer, primary_key=True, autoincrement=True, nullable=False)
    fingerprint_type_name = Column('fingerprint_type_name', String(100), nullable=False)
    
    def __repr__(self):
        return "<Fingerprint('id={}', fingerprint_type_name={} )>".format(
            self.id,
            self.fingerprint_type_name,
        )
        
class LatticeFingerprint(Base):
    __tablename__ = 'lattice_fingerprint'

    id = Column('id', Integer, primary_key=True, autoincrement=True, nullable=False)
    lattice_id = Column('lattice_id', Integer, ForeignKey("lattice.id", name="fk_lattice_lattice_id"), nullable=False)
    fingerprint_type_id = Column('fingerprint_type_id', Integer, ForeignKey("fingerprint_type.id", name="fk_fingerprint_type_fingerprint_type_id"), nullable=False)
    fingerprint = Column('fingerprint', String, nullable=False)
    lattice = relationship("Lattice", uselist=True)
    fingerprint_type = relationship("fingerprint_type", uselist=True)
     
    def __repr__(self):
        return "<LatticeFingerprint('id={}', fingerprint_type_id={}, fingerprint={} )>".format(
            self.id,
            self.fingerprint_type_id,
            self.fingerprint
        )
    