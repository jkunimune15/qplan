import pytz
import entity
import ephem

# The common solar system bodies
moon = entity.FixedTarget(name="Moon")
moon.body = ephem.Moon()
sun = entity.FixedTarget(name="Sun")
sun.body = ephem.Sun()
