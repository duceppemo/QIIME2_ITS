def format_elapsed(seconds):
    """Format a duration in seconds as a compact "1d2h3m4s"-style string."""
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
    return ''.join(f'{round(value)}{name}' for name, value in periods if value) or '0s'
