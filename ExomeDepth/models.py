#!venv/bin/python

from sqlalchemy import Column, Integer, String, UniqueConstraint

from database import Base


class Sample(Base):
    __tablename__ = 'samples'

    id = Column(Integer, primary_key=True)
    sample = Column(String, nullable=False)
    flowcell = Column(String, nullable=False)
    refset = Column(String, nullable=False)

    __table_args__ = (
        UniqueConstraint('sample', 'flowcell'),
    )

    def __repr__(self):
        return "<Sample(sample={}, flowcell={}, refset={})>".format(self.name, self.flowcell, self.refset)
