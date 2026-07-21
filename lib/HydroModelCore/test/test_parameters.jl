using HydroModelCore

@parameters p [bounds=(0, 1), description="test", guess=0.5, unit="m"]

@info getbounds(p)
@info getdescription(p)
@info getguess(p)
@info getunit(p)

isparameter(p)