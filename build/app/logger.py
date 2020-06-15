from os import environ
import logging
import inspect

class Logger():

    def __init__(self, name, level=logging.INFO):
        self.level = level if type(level) == str else getLevelName(level)
        self.log = logging.getLogger(name)
        self.log.setLevel(level if level else environ.get("LOGLEVEL", logging.INFO))

    def addConsoleHandler(self, level=None, formatted=True):
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level if level else self.level)
        if formatted and console_handler.level > 0:
            if console_handler.level >= 20:
                formatter = logging.Formatter(
                    "[ %(levelname)8s ] --- %(message)s"
                )
            else:
                formatter = logging.Formatter(
                    "%(asctime)s - %(filename)15s:%(lineno)-4s - [ %(levelname)8s ] - %(funcName)15s --- %(message)s"
                )
            console_handler.setFormatter(formatter)
        self.log.addHandler(console_handler)

    def addFileHandler(self, filename, level=None):
        file_handler = logging.FileHandler(filename)
        file_handler.setLevel(level if level else self.level)
        if file_handler.level >= 20:
            formatter = logging.Formatter(
                "[ %(levelname)8s ] --- %(message)s"
            )
        else:
            formatter = logging.Formatter(
                "%(asctime)s - %(filename)15s:%(lineno)-4s - [ %(levelname)8s ] - %(funcName)15s --- %(message)s"
            )
        file_handler.setFormatter(formatter)
        self.log.addHandler(file_handler)
    
    def getFrame(self):
        return inspect.currentframe().f_back.f_code

    def getLogger(self):
        return self.log

    def getRootLogger(self):
        return logging.getLogger()


def getLevelName(level):
    return logging.getLevelName(level)
