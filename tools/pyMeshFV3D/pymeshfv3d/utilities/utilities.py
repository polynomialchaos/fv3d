################################################################################
# @file utilities.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2021-11-23
# @copyright Copyright (c) 2021
################################################################################
import logging


class LessThanFilter(logging.Filter):
    """User defined filter for logging."""

    def __init__(self, exclusive_maximum, name=''):
        super().__init__(name=name)
        self.max_level = exclusive_maximum

    def filter(self, record):
        # non-zero return means we log this message
        return 1 if record.levelno < self.max_level else 0
