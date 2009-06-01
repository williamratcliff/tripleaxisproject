import sqlalchemy
from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey
import inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relation, backref

Base=declarative_base()

engine = create_engine('mysql://root:secret@localhost/sampleinventory')
metadata = MetaData()
meta.bind = engine

class Instrument(Base):
    __tablename__='instruments'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id=Column(Integer, primary_key=True)
    name=Column(String(50))
    def __init__(self,name):
        self.name=name
    def __repr__(self):
        return "<Instrument('%s')>"%(self.name,)

class Sample(Base):
    __tablename__='samples'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id=Column(Integer, primary_key=True)
    name=Column(String, nullable=False)
    def __init__(sample_name, chemical_formula, quantity, sample_type, hazards, current_location, 
                 comments=None, MSDS=None, date_shipped=None):
        self.sample_name=sample_name
        self.chemical_formula=chemical_formula
        self.quantity=quantity 
        self.sample_type=sample_type #(powder, crystal, liquid)
        self.hazards=hazards#(toxic, flammable, etc.)
        self.current_location=current_location
        self.comments=comments
        self.MSDS=MSDS
        self.date_shipped=date_shipped #should be a date

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

    def __init__(self,name,address1,address2,city,state,zipcode,country):
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

people_organizations=Table('person_organization',metadata,
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
    #organization_id=Column(Integer, ForeignKey('organization.id'))
    organization=relation('Organization',secondary=people_organizations, backref='people')

    
    def __init__(self,name,email,organization,phone,role='experimenter'):
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
        
        



        

        
class Date(object):
    def __init__(self,day,month, year):
        self.day=day
        self.month=month
        self.year=year

class Experiment(object):
    def __init__(date_on, date_off, flux, instrument_id, sample_id, comments=''):
        self.date_on=date_on
        self.date_off=date_off
        self.flux=flux
        self.instrument_id=instrument_id
        self.sample_id=sample_id
        
class PersonExperiment(object):
    def __init__(person_id, experiment_id):
        self.person_id=person_id
        self.experiment_id=experiment_id
        

        

        
        



        

if __name__=='__main__':
    johns_name=Name('john','smith')
    johns_address=Address()
    johns_address.address1='100 Bureau Drive'
    johns_address.address2='100 Bureau Drive'
    
    john=Person(johns_name, johns_address,'john@nist.gov')