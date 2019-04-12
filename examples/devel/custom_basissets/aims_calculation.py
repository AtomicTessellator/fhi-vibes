from hilde import Settings
from hilde.templates.aims import setup_aims
from hilde.tasks import calculate_socket

settings = Settings()

atoms = settings.atoms

calc = setup_aims(settings=settings)

calculate_socket([atoms], calc)
