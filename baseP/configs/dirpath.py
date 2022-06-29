import sys
this = sys.modules[__name__]
this.static = None

def initialize_dr(name):
    this.static = name
    static = "Locally scoped static variable. Doesn't do anything here."
