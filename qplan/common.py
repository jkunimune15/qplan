import entity
import ephem

# The common solar system bodies
moon = entity.StaticTarget(name="Moon")
moon.body = ephem.Moon()
sun = entity.StaticTarget(name="Sun")
sun.body = ephem.Sun()
