<!--
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.2
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
--->

---
---
# SOLID principles for OOP in python

https://gist.github.com/dmmeteo/f630fa04c7a79d3c132b9e9e5d037bfd
---


# S: Single responsibility
## A class should have only one job.

### <span style="color:red;">~~S~~ </span>OLID:
```
class Animal:
    def __init__(self, name: str):
        self.name = name

    def get_name(self) -> str:
        pass

    def save(self, animal: Animal):
        pass

```
The Animal class violates the SRP.
How does it violate SRP?
SRP states that classes should have one responsibility, here, we can draw out two responsibilities: animal database management and animal properties management.
The constructor and get_name manage the Animal properties while the save manages the Animal storage on a database.
How will this design cause issues in the future?
If the application changes in a way that it affects database management functions. The classes that make use of Animal properties will have to be touched and recompiled to compensate for the new changes.
You see this system smells of rigidity, it’s like a domino effect, touch one card it affects all other cards in line.
To make this conform to SRP, we create another class that will handle the sole responsibility of storing an animal to a database:


```python
class Animal:
    def __init__(self, name: str):
            self.name = name

    def get_name(self):
        pass


class AnimalDB:
    def get_animal(self) -> Animal:
        pass

    def save(self, animal: Animal):
        pass

```

# O: Open-Closed Principle
## Software entities (Classes, modules, functions) should be open for extension, not modification.

### S<span style="color:red;"> ~~O~~ </span>LID:
```
class Animal:
    def __init__(self, name: str):
        self.name = name

    def get_name(self) -> str:
        pass

animals = [
    Animal('lion'),
    Animal('mouse')
]

def animal_sound(animals: list):
    for animal in animals:
        if animal.name == 'lion':
            print('roar')

        elif animal.name == 'mouse':
            print('squeak')

animal_sound(animals)
```

The function `animal_sound` does not conform to the open-closed principle because it cannot be closed against new kinds of animals.
If we add a new animal, Snake, we have to modify the `animal_sound` function.


```python
class Animal:
    def __init__(self, name: str):
        self.name = name

    def get_name(self) -> str:
        pass

    def make_sound(self):
        pass


class Lion(Animal):

    def get_name(self) -> str:
        return self.name.title()

    def make_sound(self):
        return 'roar'


class Mouse(Animal):
    def make_sound(self):
        return 'squeak'


class Snake(Animal):
    def make_sound(self):
        return 'hiss'


def animal_sound(animals: list):
    for animal in animals:
        print(animal.make_sound())

animals = [
    Lion('lion'),
    Mouse('mouse')
]
animal_sound(animals)
```

```python
animals[0].name
animals[0].get_name()
```

# L: Liskov Substitution Principle
## A sub-class must be substitutable for its super-class.

> The aim of this principle is to ascertain that a sub-class can assume the place of its super-class without errors.
If the code finds itself checking the type of class then, it must have violated this principle.

### SO<span style="color:red;"> ~~L~~ </span>ID:
```
def animal_leg_count(animals: list):
    for animal in animals:
        if isinstance(animal, Lion):
            print(lion_leg_count(animal))
        elif isinstance(animal, Mouse):
            print(mouse_leg_count(animal))
        elif isinstance(animal, Pigeon):
            print(pigeon_leg_count(animal))

animal_leg_count(animals)
```

To make this function follow the LSP principle, we will follow this LSP requirements postulated by Steve Fenton:

* If the super-class (Animal) has a method that accepts a super-class type (Animal) parameter.
  * Its sub-class (Pigeon) should accept as argument a super-class type (Animal type) or sub-class type(Pigeon type).
* If the super-class returns a super-class type (Animal).
  * Its sub-class should return a super-class type (Animal type) or sub-class type(Pigeon).


```python
# Updated base class
class Animal:
    def __init__(self, name: str):
        self.name = name

    def get_name(self) -> str:
        pass

    def make_sound(self):
        pass

    def leg_count(self):
        pass

# Sub-classes:
class Lion(Animal):

    def get_name(self) -> str:
        return self.name.title()

    def leg_count(self):
        return 4

class Snake(Animal):
    def get_name(self) -> str:
        return self.name.title()

    def make_sound(self):
        return 'hiss'

    def leg_count(self):
        return 0


def animal_leg_count(animals: list):
    for animal in animals:
        legs = animal.leg_count() or "no"
        print(f"{animal.get_name()} the <{str(animal.__class__).split('.')[-1][:-2]}> has {legs} legs.")


animals = [
    Lion('lion'),
    Snake('snicky')
]

animal_leg_count(animals)
```

# I: ISP :: Interface Segregation Principle
## Make fine grained interfaces that are client specific.

Clients should not be forced to depend upon interfaces that they do not use.
This principle deals with the disadvantages of implementing big interfaces.

### SOL<span style="color:red;"> ~~I~~ </span>D:
Let’s look at the below IShape interface:
```
class IShape:
    def draw_square(self):
        raise NotImplementedError

    def draw_rectangle(self):
        raise NotImplementedError

    def draw_circle(self):
        raise NotImplementedError

```
This interface draws squares, circles, rectangles. class Circle, Square or Rectangle implementing the IShape
interface must define the methods draw_square(), draw_rectangle(), draw_circle().

```
class Circle(IShape):
    def draw_square(self):
        pass

    def draw_rectangle(self):
        pass

    def draw_circle(self):
        pass

class Square(IShape):
    def draw_square(self):
        pass

    def draw_rectangle(self):
        pass

    def draw_circle(self):
        pass

class Rectangle(IShape):
    def draw_square(self):
        pass

    def draw_rectangle(self):
        pass

    def draw_circle(self):
        pass

```
It’s quite funny looking at the code above. class Rectangle implements methods (draw_circle and draw_square) it has no use of,
likewise Square implementing draw_circle, and draw_rectangle, and class Circle (draw_square, draw_rectangle).

If we add another method to the IShape interface, like draw_triangle(),

```
class IShape:
    def draw_square(self):
        raise NotImplementedError

    def draw_rectangle(self):
        raise NotImplementedError

    def draw_circle(self):
        raise NotImplementedError

    def draw_triangle(self):
        raise NotImplementedError
```
the classes must implement the new method or error will be thrown.

We see that it is impossible to implement a shape that can draw a circle but not a rectangle or a square or a triangle.
We can just implement the methods to throw an error that shows the operation cannot be performed.

ISP frowns against the design of this IShape interface. clients (here Rectangle, Circle, and Square) should not
be forced to depend on methods that they do not need or use.
Also, ISP states that interfaces should perform only one job (just like the SRP principle) any extra grouping of
behavior should be abstracted away to another interface.

Here, our IShape interface performs actions that should be handled independently by other interfaces.

To make our IShape interface conform to the ISP principle, we segregate the actions to different interfaces.
Classes (Circle, Rectangle, Square, Triangle, etc) can just inherit from the IShape interface and implement their own draw behavior.

```
class IShape:
    def draw(self):
        raise NotImplementedError

class Circle(IShape):
    def draw(self):
        pass

class Square(IShape):
    def draw(self):
        pass

class Rectangle(IShape):
    def draw(self):
        pass
```
We can then use the I -interfaces to create Shape specifics like Semi Circle, Right-Angled Triangle, Equilateral Triangle, Blunt-Edged Rectangle, etc.


# D: DIP :: Dependency Inversion Principle
## Dependency should be on abstractions not concretions
  1. High-level modules should not depend upon low-level modules. Both should depend upon abstractions.
  2. Abstractions should not depend on details. Details should depend upon abstractions.

There comes a point in software development where our app will be largely composed of modules.
When this happens, we have to clear things up by using dependency injection.
High-level components depending on low-level components to function.

### SOLI<span style="color:red;"> ~~D~~ </span>:

```
class XMLHttpService(XMLHttpRequestService):
    pass

class Http:
    def __init__(self, xml_http_service: XMLHttpService):
        self.xml_http_service = xml_http_service

    def get(self, url: str, options: dict):
        self.xml_http_service.request(url, 'GET')

    def post(self, url, options: dict):
        self.xml_http_service.request(url, 'POST')
```
Here, Http is the high-level component whereas HttpService is the low-level component.
This design violates DIP A: High-level modules should not depend on low-level level modules. It should depend upon its abstraction.

This Http class is forced to depend upon the XMLHttpService class.
If we were to change to change the Http connection service, maybe we want to connect to the internet through cURL or even Mock the http service.
We will painstakingly have to move through all the instances of Http to edit the code and this violates the OCP principle.

The Http class should care less the type of Http service you are using. We make a Connection interface:
```
class Connection:
    def request(self, url: str, options: dict):
        raise NotImplementedError

```
The Connection interface has a request method. With this, we pass in an argument of type Connection to our Http class:

```
class Http:
    def __init__(self, http_connection: Connection):
        self.http_connection = http_connection

    def get(self, url: str, options: dict):
        self.http_connection.request(url, 'GET')

    def post(self, url, options: dict):
        self.http_connection.request(url, 'POST')
```

So now, no matter the type of Http connection service passed to Http it can easily connect to a network
without bothering to know the type of network connection.

We can now re-implement our XMLHttpService class to implement the Connection interface:

```
class XMLHttpService(Connection):
    xhr = XMLHttpRequest()

    def request(self, url: str, options:dict):
        self.xhr.open()
        self.xhr.send()
```

We can create many Http Connection types and pass it to our Http class without any fuss about errors.

```
class NodeHttpService(Connection):
    def request(self, url: str, options:dict):
        pass

class MockHttpService(Connection):
    def request(self, url: str, options:dict):
        pass

```

Now, we can see that both high-level modules and low-level modules depend on abstractions.
Http class(high level module) depends on the Connection interface(abstraction) and
the Http service types(low level modules) in turn, depends on the Connection interface(abstraction).

Also, this DIP will force us not to violate the Liskov Substitution Principle:
The Connection types Node-XML-MockHttpService are substitutable for their parent type Connection.

```python

```
