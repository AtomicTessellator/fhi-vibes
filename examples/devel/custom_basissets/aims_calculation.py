from vibes import Settings
from vibes.templates.aims import setup_aims
from vibes.tasks import calculate_socket

settings = Settings()

atoms = settings.atoms

calc = setup_aims(settings=settings)

calculate_socket([atoms], calc)
