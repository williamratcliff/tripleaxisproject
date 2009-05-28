import sqlalchemy
import inspect


class Name(object):
    def __init__(self,first_name, last_name, middle_initial=None, suffix=None, title=None):
        """At the moment, we imagine that these are the attributes of a name"""
        self.first_name=first_name
        self.last_name=last_name
        self.middle_initial=middle_initial
        self.suffix=suffix
        self.title=title
        
class Address(object):
    def __init__(self,address1=None,address2=None, country=None, zipcode=None, city=None, state=None):
        """This is for a US address, must be subclassed to support international addresses"""
        self.address1=address1
        self.address2=address2
        self.zipcode=zipcode
        self.city=city
        
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
        
class Instrument(object):
    def __init__(name):
        self.name=name
        
class Sample(object):
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
        
        
class Person(object):
    """A Person has a number of attributes
    
    
    """
    def __init__(self,name, address,email,organization=None, phone=None,role='experimenter'):
        """Here, The role is either experimenter, or local_contact"""
        self.name=name
        self.address=address
        self.email=email
        self.phone=phone
        self.role=role
        self.organization=organization
        
        #print inspect.getargspec(self.__init__)
        
    def __repr__(self):
        return "<Person('%s' ,'%s')>" % (self.first_name, self.last_name)


        

if __name__=='__main__':
    johns_name=Name('john','smith')
    johns_address=Address()
    johns_address.address1='100 Bureau Drive'
    johns_address.address2='100 Bureau Drive'
    
    john=Person(johns_name, johns_address,'john@nist.gov')