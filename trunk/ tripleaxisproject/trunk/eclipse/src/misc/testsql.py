import sqlalchemy
from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey, DateTime,BLOB, Text
import inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relation, backref, sessionmaker
import datetime

Base=declarative_base()

#engine = create_engine('mysql://root:secret@localhost/sampleinventory')
engine = create_engine('mysql://root:secret@localhost/sampleinventory?charset=utf8&use_unicode=0', pool_recycle=3600)
metadata = MetaData()
meta.bind = engine
Session=sessionmaker(bind=engine)
session=Session()

class Instruments(Base):
    __tablename__='instruments'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id=Column(Integer, primary_key=True)
    name=Column(String(50))
    experiments=relation(Experiments,backref=backref('instruments',
                                                     cascade='all,delete-orphan'))
    #no orphan experiments without an instrument
    def __init__(self,name):
        self.name=name
    def __repr__(self):
        return "<Instrument('%s')>"%(self.name,)

class Samples(Base):
    __tablename__='samples'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id=Column(Integer, primary_key=True)
    name=Column(String(250), nullable=False)
    chemical_formula=Column(String(50))
    quantity=Column(String(50))
    sample_type=Column(String(50)) #powder, crystal, liquid, etc.
    hazards=Column(String(250))  #flamable, toxic, etc
    MSDS=Column(Text)
    date_shipped=Column(DateTime)
    current_location=Column(Text)
    comments=Column(Text)
    previous_id=Column(Integer)
    
    photos=relation(Photos,backref=backref('samples',
                                           cascade='all,delete-orphan'))
    instruments=relation(Instruments,backref=backref('samples'))
    experiments=relation(Experiments,backref=backref('samples',
                                                     cascade='all,delete-orphan'))
  
    # no photos, or experiments without a sample
    def __init__(self,sample_name, chemical_formula, quantity, sample_type, hazards, current_location, 
                 comments=None, MSDS=None, date_shipped=None):
        self.name=sample_name #eventually a barcode
        self.chemical_formula=chemical_formula
        self.quantity=quantity 
        self.sample_type=sample_type #(powder, crystal, liquid)
        self.hazards=hazards#(toxic, flammable, etc.)
        self.current_location=current_location
        self.comments=comments
        self.MSDS=MSDS
        self.date_shipped=date_shipped #should be a date

class Photos(Base):
    __tablename__='photos'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id=Column(Integer, primary_key=True)
    name=Column(String, nullable=False)
    sample_id=Column(Integer, ForeignKey('samples.id'))
    sample=relation(Samples,backref=backref('photos'))#,
                                            #cascade='all,delete-orphan'))
    experiment_id=Column(Integer, ForeignKey('experiments.id'))
    experiment=relation(Experiments,backref=backref('photos'))#,
                        #cascade='all,delete-orphan'))
    #experiment=relation('Experiments',secondary=photos_experiments, backref='photos')
    comments=Column(Text)
    def __init__(self,photo):
        self.photo=photo

#photos_experiments=Table('photos_experiment',metadata,
#                          Column('photo_id',Integer,ForeignKey('photos.id')),
#                          Column('experiment_id',Integer,ForeignKey('experiments.id'))
#                          )        


class Experiments(object):
    __tablename__='experiments'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id=Column(Integer, primary_key=True)
    name=Column(String(250), nullable=False) #Title from Webpage
    proposal_id=Column(Integer,nullable=False)
    date_on=Column(DateTime)
    date_off=Column(DateTime)
    flux=Column(Integer)
    comments=Column(Text)
    sample_id=Column(Integer, ForeignKey('samples.id'))
    sample=relation(Samples,backref=backref('experiments'))
    instrument_id=Column(Integer, ForeignKey('experiments.id'))
    instrument=relation(Instruments,backref=backref('experiments'))
    
    
    def __init__(date_on, date_off, flux, instrument_id, sample_id, comments=''):
        self.date_on=date_on
        self.date_off=date_off
        self.flux=flux
        self.instrument_id=instrument_id
        self.sample_id=sample_id

people_experiments=Table('people_experiments',metadata,
                          Column('person_id',Integer,ForeignKey('people.id')),
                          Column('experiment_id',Integer,ForeignKey('experiments.id'))
                          )        

class Address(object):
    def __init__(self,address1=None,address2=None, country=None, zipcode=None, city=None, state=None):
        """This is for a US address, must be subclassed to support international addresses"""
        self.address1=address1
        self.address2=address2
        self.zipcode=zipcode
        self.city=city
        self.state=state
        self.country=country
        
class Organization(Base):
    """An Organization has a number of attributes
    """
    __tablename__='organizations'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id=Column(Integer, primary_key=True)
    name=Column(String,nullable=False)
    address1=Column(String, nullable=False)
    address2=Column(String)
    city=Column(String)
    state=Column(String)
    zipcode=Column(Integer)
    country=Column(String, nullable=False)        

    def __init__(self,name,address1=None,address2=None,city=None,state=None,zipcode=None,country='USA'):
        """The properties of an organization"""
        self.name=name
        self.address1=address1
        self.address2=address2
        self.city=city
        self.state=state
        self.zipcode=zipcode
        self.country=country
        
class Name(object):
    def __init__(self,first_name, last_name, middle_initial=None, suffix=None, title=None):
        """At the moment, we imagine that these are the attributes of a name"""
        self.first_name=first_name
        self.last_name=last_name
        self.middle_initial=middle_initial
        self.suffix=suffix
        self.title=title

people_organizations=Table('people_organization',metadata,
                          Column('person_id',Integer,ForeignKey('people.id')),
                          Column('organization_id',Integer,ForeignKey('organizations.id'))
                          )
        
class Person(Base):
    """A Person has a number of attributes
    """
    __tablename__='people'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id=Column(Integer, primary_key=True)
    title=Column(String)
    first_name=Column(String, nullable=False)
    middle_initial=Column(String)
    last_name=Column(String, nullable=False)
    suffix=Column(String)
    email=Column(String, nullable=False)
    phone=Column(String, nullable=False)
    role=Column(String, nullable=False)
    ##organization_id=Column(Integer, ForeignKey('organization.id'))
    organization=relation('Organization',secondary=people_organizations, backref='people')

    
    def __init__(self,name,email=None,organization=None,phone=None,role='experimenter'):
        """Here, The role is either experimenter, or local_contact"""
        self.first_name=name.first_name
        self.middle_initial=name.middle_initial
        self.last_name=name.last_name
        self.title=name.title
        self.suffix=name.suffix
        self.email=email
        self.phone=phone
        self.role=role
        self.organization=organization
        
        #print inspect.getargspec(self.__init__)
        
    def __repr__(self):
        return "<Person('%s' ,'%s')>" % (self.first_name, self.last_name)        
        
        

        

if __name__=='__main__':
    williams_name=Name('william','ratcliff')
    NCNR=Organization('NCNR')
    NCNR.address1='100 Bureau Drive'
    NCNR.address2='MS 610, Rm E151'
    NCNR.city='Gaithersburg'
    NCNR.state='md'
    NCNR.zipcode=20886
    william=Person(williams_name)
    william.email='william.ratcliff@nist.gov'
    william.phone='301-975-4316'
    william.role='nist_staff'
    william.organization=NCNR
    