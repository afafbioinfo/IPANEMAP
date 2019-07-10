import time
import sys
import os

class Progress:
    def __init__(self):
        self.tasks = []
        self.immediate = False
        self.limitTitleDepth = 3

    def Print(self, txt, endl = True, output=False):
        if self.immediate:
            sys.stderr.write("]"+os.linesep)
        sys.stderr.write(self.getIndentTxt())
        if output:
            sys.stdout.write(txt)
            if endl:
                sys.stdout.write(os.linesep)
        else:
            sys.stderr.write(txt)
            if endl:
                sys.stderr.write(os.linesep)
        self.immediate = False


    def StartTask(self, title):
        self.Print("[%s..."%title, False)
        self.tasks.append((title, time.time()))
        self.immediate = True


    def EndTask(self, showTime = True):
        (t, tim) = self.tasks.pop()
        added = ""
        elapsedTime = time.time()-tim
        if showTime and elapsedTime>1.:
            added = " in %.2f secs"%(elapsedTime)
        if self.immediate:
            sys.stderr.write(" Done%s]"%(added)+os.linesep)
        else:
            sys.stderr.write(self.getIndentTxt() + "[Done%s]"%(added) + os.linesep)
        self.immediate = False

    def Flush(self):
        while len(self.tasks)>0:
            self.EndTask()

    def getIndentTxt(self):
        return "  "*len(self.tasks)

progress = Progress()
