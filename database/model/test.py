from sqlalchemy import Column, Integer, String, Float, DateTime
from ..settings import Engine
from ..settings import Base

class Test(Base):
    """
    テストモデル
    """

    __tablename__ = 'tests'
    __table_args__ = {
        'comment': 'DB接続テスト'
    }

    id = Column('id', Integer, primary_key=True, autoincrement=True)
    content = Column('content', String(200))
    title = Column('title', String(100))

if __name__ == "__main__":
    Base.metadata.create_all(bind=Engine)
