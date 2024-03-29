"""Notifications via, e.g., email"""


def send_simple_mail(message: str, to_addr: str, extra_message: str = None) -> None:
    """
    Send simple e-mail message

    Args:
    ----
      message: the message
      to_addr: the mail address of recipient
      extra_message:  prefix of message

    """
    import os

    log = os.system(
        f'echo "{extra_message}" | mailx -s "[vibes] {message:s}" {to_addr:s}'
    )

    if log:
        print(f"Sending the Mail returned error code {log!s:s}")
