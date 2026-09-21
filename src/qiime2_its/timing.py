def format_elapsed(seconds):
    """Format a duration in seconds as a compact "1d2h3m4s"-style string."""
    # Rounded once, up front: rounding each unit after the split turned
    # 119.7 s into "1m60s".
    minutes, seconds = divmod(round(seconds), 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
    return ''.join(f'{value}{name}' for name, value in periods if value) or '0s'
