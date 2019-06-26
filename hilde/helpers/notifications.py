""" Notifications via, e.g., email """


def send_simple_mail(message, to_addr, extra_message=""):
    """Send simple e-mail message

    Parameters
    ----------
    message: str
        Message to send
    to_addr: str
        The e-mail address
    extra_message:
        Additional message to send
    """
    import os

    log = os.system(
        'echo "{}" | mailx -s "[hilde] {:s}" {:s}'.format(
            extra_message, message, to_addr
        )
    )

    if log:
        print("Sending the Mail returned error code {:s}".format(str(log)))
