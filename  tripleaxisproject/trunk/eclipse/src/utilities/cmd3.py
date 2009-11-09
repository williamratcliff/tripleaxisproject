"""Variant on standard library's cmd with extra features.

To use, simply import cmd2.Cmd instead of cmd.Cmd; use precisely as though you
were using the standard library's cmd, while enjoying the extra features.

Searchable command history (commands: "hi", "li", "run")
Load commands from file, save to file, edit commands in file
Multi-line commands
Case-insensitive commands
Special-character shortcut commands (beyond cmd's "@" and "!")
Settable environment parameters
Optional _onchange_{paramname} called when environment parameter changes
Parsing commands with `optparse` options (flags)
Redirection to file with >, >>; input from file with <
Easy transcript-based testing of applications (see example/example.py)

Note that redirection with > and | will only work if `self.stdout.write()`
is used in place of `print`.  The standard library's `cmd` module is 
written to use `self.stdout.write()`, 

- Catherine Devlin, Jan 03 2008 - catherinedevlin.blogspot.com

mercurial repository at http://www.assembla.com/wiki/show/python-cmd2
"""
import cmd, re, os, sys, optparse, subprocess, tempfile, pyparsing, doctest
import unittest, string, datetime, urllib, glob
from code import InteractiveConsole, InteractiveInterpreter, softspace
from optparse import make_option
import sequenceloader
__version__ = '0.5.5'

class OptionParser(optparse.OptionParser):
    def exit(self, status=0, msg=None):
        self.values._exit = True
        if msg:
            print msg
            
    def print_help(self, *args, **kwargs):
        # now, I need to call help of the calling function.  Hmm.
        try:
            print self._func.__doc__
        except AttributeError:
            pass
        optparse.OptionParser.print_help(self, *args, **kwargs)

    def error(self, msg):
        """error(msg : string)

        Print a usage message incorporating 'msg' to stderr and exit.
        If you override this in a subclass, it should not return -- it
        should either exit or raise an exception.
        """
        raise
        
def remainingArgs(oldArgs, newArgList):
    '''
    >>> remainingArgs('-f bar   bar   cow', ['bar', 'cow'])
    'bar   cow'
    '''
    pattern = '\s+'.join(re.escape(a) for a in newArgList) + '\s*$'
    matchObj = re.search(pattern, oldArgs)
    return oldArgs[matchObj.start():]

def options(option_list):
    def option_setup(func):
        optionParser = OptionParser()
        for opt in option_list:
            optionParser.add_option(opt)
        optionParser.set_usage("%s [options] arg" % func.__name__.strip('do_'))
        optionParser._func = func
        def newFunc(instance, arg):
            try:
                opts, newArgList = optionParser.parse_args(arg.split()) # doesn't understand quoted strings shouldn't be dissected!
                newArgs = remainingArgs(arg, newArgList)  # should it permit flags after args?
            except (optparse.OptionValueError, optparse.BadOptionError,
                    optparse.OptionError, optparse.AmbiguousOptionError,
                    optparse.OptionConflictError), e:
                print e
                optionParser.print_help()
                return
            if hasattr(opts, '_exit'):
                return None
            if hasattr(arg, 'parser'):
                terminator = arg.parsed.terminator
                try:
                    if arg.parsed.terminator[0] == '\n':
                        terminator = arg.parsed.terminator[0]
                except IndexError:
                    pass
                arg = arg.parser('%s %s%s%s' % (arg.parsed.command, newArgs, 
                                                terminator, arg.parsed.suffix))
            else:
                arg = newArgs
            result = func(instance, arg, opts)                            
            return result        
        newFunc.__doc__ = '%s\n%s' % (func.__doc__, optionParser.format_help())
        return newFunc
    return option_setup

class PasteBufferError(EnvironmentError):
    if sys.platform[:3] == 'win':
        errmsg = """Redirecting to or from paste buffer requires pywin32
to be installed on operating system.
Download from http://sourceforge.net/projects/pywin32/"""
    else:
        errmsg = """Redirecting to or from paste buffer requires xclip 
to be installed on operating system.
On Debian/Ubuntu, 'sudo apt-get install xclip' will install it."""        
    def __init__(self):
        Exception.__init__(self, self.errmsg)

'''check here if functions exist; otherwise, stub out'''
pastebufferr = """Redirecting to or from paste buffer requires %s
to be installed on operating system.
%s"""
if subprocess.mswindows:
    try:
        import win32clipboard
        def getPasteBuffer():
            win32clipboard.OpenClipboard(0)
            try:
                result = win32clipboard.GetClipboardData()
            except TypeError:
                result = ''  #non-text
            win32clipboard.CloseClipboard()
            return result            
        def writeToPasteBuffer(txt):
            win32clipboard.OpenClipboard(0)
            win32clipboard.EmptyClipboard()
            win32clipboard.SetClipboardText(txt)
            win32clipboard.CloseClipboard()        
    except ImportError:
        def getPasteBuffer(*args):
            raise OSError, pastebufferr % ('pywin32', 'Download from http://sourceforge.net/projects/pywin32/')
        setPasteBuffer = getPasteBuffer
else:
    can_clip = False
    try:
        subprocess.check_call('xclip -o -sel clip', shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        can_clip = True
    except AttributeError:  # check_call not defined, Python < 2.5
        teststring = 'Testing for presence of xclip.'
        xclipproc = subprocess.Popen('xclip -sel clip', shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        xclipproc.stdin.write(teststring)
        xclipproc.stdin.close()
        xclipproc = subprocess.Popen('xclip -o -sel clip', shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)        
        if xclipproc.stdout.read() == teststring:
            can_clip = True
    except (subprocess.CalledProcessError, OSError, IOError):
        pass
    if can_clip:    
        def getPasteBuffer():
            xclipproc = subprocess.Popen('xclip -o -sel clip', shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
            return xclipproc.stdout.read()
        def writeToPasteBuffer(txt):
            xclipproc = subprocess.Popen('xclip -sel clip', shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
            xclipproc.stdin.write(txt)
            xclipproc.stdin.close()
            # but we want it in both the "primary" and "mouse" clipboards
            xclipproc = subprocess.Popen('xclip', shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
            xclipproc.stdin.write(txt)
            xclipproc.stdin.close()
    else:
        def getPasteBuffer(*args):
            raise OSError, pastebufferr % ('xclip', 'On Debian/Ubuntu, install with "sudo apt-get install xclip"')
        setPasteBuffer = getPasteBuffer
        writeToPasteBuffer = getPasteBuffer
          
pyparsing.ParserElement.setDefaultWhitespaceChars(' \t')

class ParsedString(str):
    pass

class SkipToLast(pyparsing.SkipTo):
    def parseImpl( self, instring, loc, doActions=True ):
        resultStore = []
        startLoc = loc
        instrlen = len(instring)
        expr = self.expr
        failParse = False
        while loc <= instrlen:
            try:
                if self.failOn:
                    failParse = True
                    self.failOn.tryParse(instring, loc)
                    failParse = False
                loc = expr._skipIgnorables( instring, loc )
                expr._parse( instring, loc, doActions=False, callPreParse=False )
                skipText = instring[startLoc:loc]
                if self.includeMatch:
                    loc,mat = expr._parse(instring,loc,doActions,callPreParse=False)
                    if mat:
                        skipRes = ParseResults( skipText )
                        skipRes += mat
                        resultStore.append((loc, [ skipRes ]))
                    else:
                        resultStore,append((loc, [ skipText ]))
                else:
                    resultStore.append((loc, [ skipText ]))
                loc += 1
            except (pyparsing.ParseException,IndexError):
                if failParse:
                    raise
                else:
                    loc += 1
        if resultStore:
            return resultStore[-1]
        else:
            exc = self.myException
            exc.loc = loc
            exc.pstr = instring
            raise exc    

def replace_with_file_contents(fname):
    if fname:
        try:
            result = open(os.path.expanduser(fname[0])).read()
        except IOError:
            result = '< %s' % fname[0]  # wasn't a file after all
    else:
        result = getPasteBuffer()
    return result      

class EmbeddedConsoleExit(Exception):
    pass

class MyInteractiveConsole(InteractiveConsole):
    def runcode(self, code):
        """Execute a code object.

        When an exception occurs, self.showtraceback() is called to
        display a traceback.  All exceptions are caught except
        SystemExit, which is reraised.

        A note about KeyboardInterrupt: this exception may occur
        elsewhere in this code, and may not always be caught.  The
        caller should be prepared to deal with it.
        
        Copied directly from code.InteractiveInterpreter, except for
        EmbeddedConsoleExit exceptions.
        """
        try:
            exec code in self.locals
        except (SystemExit, EmbeddedConsoleExit):
            raise
        except:
            self.showtraceback()
        else:
            if softspace(sys.stdout, 0):
                print

class Cmd(cmd.Cmd):
    echo = False
    case_insensitive = True
    continuation_prompt = '> '  
    timing = False
    legalChars = '!#$%.:?@_' + pyparsing.alphanums + pyparsing.alphas8bit  # make sure your terminators are not in here!
    shortcuts = {'?': 'help', '!': 'shell', '@': 'load', '@@': '_relative_load'}
    excludeFromHistory = '''run r list l history hi ed edit li eof'''.split()
    noSpecialParse = 'set ed edit exit'.split()
    defaultExtension = 'txt'
    default_file_name = 'command.txt'
    abbrev = True
    nonpythoncommand = 'cmd'
    current_script_dir = None
    reserved_words = []
    feedback_to_output = False
    quiet = False
    settable = ['prompt', 'continuation_prompt', 'default_file_name', 'editor',
                'case_insensitive', 'feedback_to_output', 'quiet', 'echo', 'timing', 
                'abbrev']   
    settable.sort()
    
    def poutput(self, msg):
        self.stdout.write(msg)
        if msg[-1] != '\n':
            self.stdout.write('\n')
    def perror(self, errmsg, statement=None):
        print str(errmsg)
    def pfeedback(self, msg):
        """For printing nonessential feedback.  Can be silenced with `quiet`.
           Inclusion in redirected output is controlled by `feedback_to_output`."""
        if not self.quiet:
            if self.feedback_to_output:
                self.poutput(msg)
            else:
                print msg
    editor = os.environ.get('EDITOR')
    _STOP_AND_EXIT = 2
    if not editor:
        if sys.platform[:3] == 'win':
            editor = 'notepad'
        else:
            for editor in ['gedit', 'kate', 'vim', 'emacs', 'nano', 'pico']:
                if not os.system('which %s' % (editor)):
                    break
            
    def do_cmdenvironment(self, args):
        '''Summary report of interactive parameters.'''
        self.stdout.write("""
        Commands are %(casesensitive)scase-sensitive.
        Commands may be terminated with: %(terminators)s
        Settable parameters: %(settable)s\n""" % \
        { 'casesensitive': (self.case_insensitive and 'not ') or '',
          'terminators': str(self.terminators),
          'settable': ' '.join(self.settable)
        })
        
    def do_help(self, arg):
        if arg:
            funcname = self.func_named(arg)
            if funcname:
                fn = getattr(self, funcname)
                try:
                    fn.optionParser.print_help(file=self.stdout)
                except AttributeError:
                    cmd.Cmd.do_help(self, funcname[3:])
        else:
            cmd.Cmd.do_help(self, arg)
        
    def __init__(self, *args, **kwargs):        
        cmd.Cmd.__init__(self, *args, **kwargs)
        self.initial_stdout = sys.stdout
        self.history = History()
        self._init_parser()
        self.pystate = {}
        self.shortcuts = sorted(self.shortcuts.items(), reverse=True)
        self.keywords = self.reserved_words + [fname[3:] for fname in dir(self) 
                                               if fname.startswith('do_')]        
    def do_shortcuts(self, args):
        """Lists single-key shortcuts available."""
        result = "\n".join('%s: %s' % (sc[0], sc[1]) for sc in sorted(self.shortcuts))
        self.stdout.write("Single-key shortcuts for other commands:\n%s\n" % (result))

    prefixParser = pyparsing.Empty()
    commentGrammars = pyparsing.Or([pyparsing.pythonStyleComment, pyparsing.cStyleComment])
    commentGrammars.addParseAction(lambda x: '')
    commentInProgress  = pyparsing.Literal('/*') + pyparsing.SkipTo(pyparsing.stringEnd)
    commentInProgress = pyparsing.NoMatch()
    terminators = [';']
    blankLinesAllowed = False
    multilineCommands = []
    
    def _init_parser(self):
        r'''
        >>> c = Cmd()
        >>> c.multilineCommands = ['multiline']
        >>> c.case_insensitive = True
        >>> c._init_parser()
        >>> print c.parser.parseString('').dump()        
        []        
        >>> print c.parser.parseString('/* empty command */').dump()        
        []        
        >>> print c.parser.parseString('plainword').dump()        
        ['plainword', '']
        - command: plainword
        - statement: ['plainword', '']
          - command: plainword        
        >>> print c.parser.parseString('termbare;').dump()
        ['termbare', '', ';', '']
        - command: termbare
        - statement: ['termbare', '', ';']
          - command: termbare
          - terminator: ;
        - terminator: ;        
        >>> print c.parser.parseString('termbare; suffx').dump()
        ['termbare', '', ';', 'suffx']
        - command: termbare
        - statement: ['termbare', '', ';']
          - command: termbare
          - terminator: ;
        - suffix: suffx
        - terminator: ;        
        >>> print c.parser.parseString('barecommand').dump()
        ['barecommand', '']
        - command: barecommand
        - statement: ['barecommand', '']
          - command: barecommand
        >>> print c.parser.parseString('COMmand with args').dump()
        ['command', 'with args']
        - args: with args
        - command: command
        - statement: ['command', 'with args']
          - args: with args
          - command: command
        >>> print c.parser.parseString('command with args and terminator; and suffix').dump()
        ['command', 'with args and terminator', ';', 'and suffix']
        - args: with args and terminator
        - command: command
        - statement: ['command', 'with args and terminator', ';']
          - args: with args and terminator
          - command: command
          - terminator: ;
        - suffix: and suffix
        - terminator: ;
        >>> print c.parser.parseString('simple | piped').dump()
        ['simple', '', '|', ' piped']
        - command: simple
        - pipeTo:  piped
        - statement: ['simple', '']
          - command: simple
        >>> print c.parser.parseString('double-pipe || is not a pipe').dump()
        ['double', '-pipe || is not a pipe']
        - args: -pipe || is not a pipe
        - command: double
        - statement: ['double', '-pipe || is not a pipe']
          - args: -pipe || is not a pipe
          - command: double
        >>> print c.parser.parseString('command with args, terminator;sufx | piped').dump()
        ['command', 'with args, terminator', ';', 'sufx', '|', ' piped']
        - args: with args, terminator
        - command: command
        - pipeTo:  piped
        - statement: ['command', 'with args, terminator', ';']
          - args: with args, terminator
          - command: command
          - terminator: ;
        - suffix: sufx
        - terminator: ;
        >>> print c.parser.parseString('output into > afile.txt').dump()
        ['output', 'into', '>', 'afile.txt']
        - args: into
        - command: output
        - output: >
        - outputTo: afile.txt
        - statement: ['output', 'into']
          - args: into
          - command: output   
        >>> print c.parser.parseString('output into;sufx | pipethrume plz > afile.txt').dump()
        ['output', 'into', ';', 'sufx', '|', ' pipethrume plz', '>', 'afile.txt']
        - args: into
        - command: output
        - output: >
        - outputTo: afile.txt
        - pipeTo:  pipethrume plz
        - statement: ['output', 'into', ';']
          - args: into
          - command: output
          - terminator: ;
        - suffix: sufx
        - terminator: ;
        >>> print c.parser.parseString('output to paste buffer >> ').dump()
        ['output', 'to paste buffer', '>>', '']
        - args: to paste buffer
        - command: output
        - output: >>
        - statement: ['output', 'to paste buffer']
          - args: to paste buffer
          - command: output
        >>> print c.parser.parseString('ignore the /* commented | > */ stuff;').dump()
        ['ignore', 'the /* commented | > */ stuff', ';', '']
        - args: the /* commented | > */ stuff
        - command: ignore
        - statement: ['ignore', 'the /* commented | > */ stuff', ';']
          - args: the /* commented | > */ stuff
          - command: ignore
          - terminator: ;
        - terminator: ;
        >>> print c.parser.parseString('has > inside;').dump()
        ['has', '> inside', ';', '']
        - args: > inside
        - command: has
        - statement: ['has', '> inside', ';']
          - args: > inside
          - command: has
          - terminator: ;
        - terminator: ;        
        >>> print c.parser.parseString('multiline has > inside an unfinished command').dump()      
        ['multiline', ' has > inside an unfinished command']
        - multilineCommand: multiline        
        >>> print c.parser.parseString('multiline has > inside;').dump()
        ['multiline', 'has > inside', ';', '']
        - args: has > inside
        - multilineCommand: multiline
        - statement: ['multiline', 'has > inside', ';']
          - args: has > inside
          - multilineCommand: multiline
          - terminator: ;
        - terminator: ;        
        >>> print c.parser.parseString('multiline command /* with comment in progress;').dump()
        ['multiline', ' command']
        - multilineCommand: multiline
        >>> print c.parser.parseString('multiline command /* with comment complete */ is done;').dump()
        ['multiline', 'command /* with comment complete */ is done', ';', '']
        - args: command /* with comment complete */ is done
        - multilineCommand: multiline
        - statement: ['multiline', 'command /* with comment complete */ is done', ';']
          - args: command /* with comment complete */ is done
          - multilineCommand: multiline
          - terminator: ;
        - terminator: ;
        >>> print c.parser.parseString('multiline command ends\n\n').dump()
        ['multiline', 'command ends', '\n', '\n']
        - args: command ends
        - multilineCommand: multiline
        - statement: ['multiline', 'command ends', '\n', '\n']
          - args: command ends
          - multilineCommand: multiline
          - terminator: ['\n', '\n']
        - terminator: ['\n', '\n']
        '''
        outputParser = (pyparsing.Literal('>>') | (pyparsing.WordStart() + '>') | pyparsing.Regex('[^=]>'))('output')
        
        terminatorParser = pyparsing.Or([(hasattr(t, 'parseString') and t) or pyparsing.Literal(t) for t in self.terminators])('terminator')
        stringEnd = pyparsing.stringEnd ^ '\nEOF'
        self.multilineCommand = pyparsing.Or([pyparsing.Keyword(c, caseless=self.case_insensitive) for c in self.multilineCommands])('multilineCommand')
        oneLineCommand = (~self.multilineCommand + pyparsing.Word(self.legalChars))('command')
        pipe = pyparsing.Keyword('|', identChars='|')
        self.commentGrammars.ignore(pyparsing.quotedString).setParseAction(lambda x: '')
        self.commentInProgress.ignore(pyparsing.quotedString).ignore(pyparsing.cStyleComment)       
        afterElements = \
            pyparsing.Optional(pipe + pyparsing.SkipTo(outputParser ^ stringEnd)('pipeTo')) + \
            pyparsing.Optional(outputParser + pyparsing.SkipTo(stringEnd).setParseAction(lambda x: x[0].strip())('outputTo'))
        if self.case_insensitive:
            self.multilineCommand.setParseAction(lambda x: x[0].lower())
            oneLineCommand.setParseAction(lambda x: x[0].lower())
        if self.blankLinesAllowed:
            self.blankLineTerminationParser = pyparsing.NoMatch
        else:
            self.blankLineTerminator = (pyparsing.lineEnd + pyparsing.lineEnd)('terminator')
            self.blankLineTerminator.setResultsName('terminator')
            self.blankLineTerminationParser = ((self.multilineCommand ^ oneLineCommand) + pyparsing.SkipTo(self.blankLineTerminator).setParseAction(lambda x: x[0].strip())('args') + self.blankLineTerminator)('statement')
        self.multilineParser = (((self.multilineCommand ^ oneLineCommand) + SkipToLast(terminatorParser).setParseAction(lambda x: x[0].strip())('args') + terminatorParser)('statement') +
                                pyparsing.SkipTo(outputParser ^ pipe ^ stringEnd).setParseAction(lambda x: x[0].strip())('suffix') + afterElements)
        self.multilineParser.ignore(self.commentInProgress)
        self.singleLineParser = ((oneLineCommand + pyparsing.SkipTo(terminatorParser ^ stringEnd ^ pipe ^ outputParser).setParseAction(lambda x:x[0].strip())('args'))('statement') +
                                 pyparsing.Optional(terminatorParser) + afterElements)
        self.multilineParser = self.multilineParser.setResultsName('multilineParser')
        self.singleLineParser = self.singleLineParser.setResultsName('singleLineParser')
        self.blankLineTerminationParser = self.blankLineTerminationParser.setResultsName('blankLineTerminatorParser')
        self.parser = self.prefixParser + (
            stringEnd |
            self.multilineParser |
            self.singleLineParser |
            self.blankLineTerminationParser | 
            self.multilineCommand + pyparsing.SkipTo(stringEnd)
            )
        self.parser.ignore(pyparsing.quotedString).ignore(self.commentGrammars)
        
        inputMark = pyparsing.Literal('<')
        inputMark.setParseAction(lambda x: '')
        fileName = pyparsing.Word(self.legalChars + '/\\')
        inputFrom = fileName('inputFrom')
        inputFrom.setParseAction(replace_with_file_contents)
        # a not-entirely-satisfactory way of distinguishing < as in "import from" from <
        # as in "lesser than"
        self.inputParser = inputMark + pyparsing.Optional(inputFrom) + pyparsing.Optional('>') + \
                           pyparsing.Optional(fileName) + (pyparsing.stringEnd | '|')
        self.inputParser.ignore(pyparsing.quotedString).ignore(self.commentGrammars).ignore(self.commentInProgress)               
    
    def preparse(self, raw, **kwargs):
        return raw
    
    def parsed(self, raw, **kwargs):
        if isinstance(raw, ParsedString):
            p = raw
        else:
            raw = self.preparse(raw, **kwargs)
            s = self.inputParser.transformString(raw.lstrip())
            for (shortcut, expansion) in self.shortcuts:
                if s.lower().startswith(shortcut):
                    s = s.replace(shortcut, expansion + ' ', 1)
                    break
            result = self.parser.parseString(s)
            result['command'] = result.multilineCommand or result.command        
            result['raw'] = raw
            result['clean'] = self.commentGrammars.transformString(result.args)
            result['expanded'] = s        
            p = ParsedString(result.clean)
            p.parsed = result
            p.parser = self.parsed
        for (key, val) in kwargs.items():
            p.parsed[key] = val
        return p
              
    def postparsing_precmd(self, statement):
        stop = 0
        return stop, statement

    def postparsing_postcmd(self, stop):
        return stop
    def func_named(self, arg):
        result = None
        target = 'do_' + arg
        if target in dir(self):
            result = target
        else:
            if self.abbrev:   # accept shortened versions of commands
                funcs = [fname for fname in self.keywords if fname.startswith(arg)]
                if len(funcs) == 1:
                    result = 'do_' + funcs[0]
        return result
    def onecmd(self, line):
        """Interpret the argument as though it had been typed in response
        to the prompt.

        This may be overridden, but should not normally need to be;
        see the precmd() and postcmd() methods for useful execution hooks.
        The return value is a flag indicating whether interpretation of
        commands by the interpreter should stop.
        
        This (`cmd2`) version of `onecmd` already override's `cmd`'s `onecmd`.

        """
        if not line:
            return self.emptyline()
        if not pyparsing.Or(self.commentGrammars).setParseAction(lambda x: '').transformString(line):
            return 0    # command was empty except for comments
        try:
            statement = self.parsed(line)
            while statement.parsed.multilineCommand and (statement.parsed.terminator == ''):
                statement = '%s\n%s' % (statement.parsed.raw, 
                                        self.pseudo_raw_input(self.continuation_prompt))                
                statement = self.parsed(statement)
        except Exception, e:
            self.perror(e)
            return 0
        if statement.parsed.command not in self.excludeFromHistory:
            self.history.append(statement.parsed.raw)
        try:
            (stop, statement) = self.postparsing_precmd(statement)
        except Exception, e:
            self.perror(e)
            return 0
        if stop:
            return self.postparsing_postcmd(stop)
        
        if not statement.parsed.command:
            return self.postparsing_postcmd(stop=0)
                
        statekeeper = None

        if statement.parsed.pipeTo:
            redirect = subprocess.Popen(statement.parsed.pipeTo, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
            statekeeper = Statekeeper(self, ('stdout',))   
            self.stdout = redirect.stdin
        elif statement.parsed.output:
            statekeeper = Statekeeper(self, ('stdout',))            
            if statement.parsed.outputTo:
                mode = 'w'
                if statement.parsed.output == '>>':
                    mode = 'a'
                try:
                    self.stdout = open(os.path.expanduser(statement.parsed.outputTo), mode)                            
                except OSError, e:
                    self.perror(e)
                    return self.postparsing_postcmd(stop=0) 
            else:
                statekeeper = Statekeeper(self, ('stdout',))
                self.stdout = tempfile.TemporaryFile()
                if statement.parsed.output == '>>':
                    self.stdout.write(getPasteBuffer())
        try:
            try:
                # "heart" of the command, replace's cmd's onecmd()
                self.lastcmd = statement.parsed.expanded   
                funcname = self.func_named(statement.parsed.command)
                if not funcname:
                    return self.postparsing_postcmd(self.default(statement))  
                try:
                    func = getattr(self, funcname)
                except AttributeError:
                    return self.postparsing_postcmd(self.default(statement))                  
                timestart = datetime.datetime.now()
                stop = func(statement) 
                if self.timing:
                    self.pfeedback('Elapsed: %s' % str(datetime.datetime.now() - timestart))
            except Exception, e:
                self.perror(e, statement)
        finally:
            if statekeeper:
                if statement.parsed.output and not statement.parsed.outputTo:
                    self.stdout.seek(0)
                    try:
                        writeToPasteBuffer(self.stdout.read())
                    except Exception, e:
                        self.perror(e)
                elif statement.parsed.pipeTo:
                    for result in redirect.communicate():              
                        statekeeper.stdout.write(result or '')                        
                self.stdout.close()
                statekeeper.restore()
                                 
            return self.postparsing_postcmd(stop)        
        
    def pseudo_raw_input(self, prompt):
        """copied from cmd's cmdloop; like raw_input, but accounts for changed stdin, stdout"""
        
        if self.use_rawinput:
            try:
                line = raw_input(prompt)
            except EOFError:
                line = 'EOF'
        else:
            self.stdout.write(prompt)
            self.stdout.flush()
            line = self.stdin.readline()
            if not len(line):
                line = 'EOF'
            else:
                if line[-1] == '\n': # this was always true in Cmd
                    line = line[:-1] 
        return line
                          
    def cmdloop(self, intro=None):
        """Repeatedly issue a prompt, accept input, parse an initial prefix
        off the received input, and dispatch to action methods, passing them
        the remainder of the line as argument.
        """

        # An almost perfect copy from Cmd; however, the pseudo_raw_input portion
        # has been split out so that it can be called separately
        
        sys.stderr.write('doing something\n')
        sys.stderr.flush()
        self.preloop()
        if self.use_rawinput and self.completekey:
            try:
                import readline
                self.old_completer = readline.get_completer()
                readline.set_completer(self.complete)
                readline.parse_and_bind(self.completekey+": complete")
            except ImportError:
                sys.stdout.write('cannot import readline\n')
                sys.stdout.flush()
                pass
        try:
            if intro is not None:
                self.intro = intro
            if self.intro:
                self.stdout.write(str(self.intro)+"\n")
            stop = None
            while not stop:
                if self.cmdqueue:
                    line = self.cmdqueue.pop(0)
                else:
                    line = self.pseudo_raw_input(self.prompt)
                if (self.echo) and (isinstance(self.stdin, file)):
                    self.stdout.write(line + '\n')
                line = self.precmd(line)
                stop = self.onecmd(line)
                stop = self.postcmd(stop, line)
            self.postloop()
        finally:
            if self.use_rawinput and self.completekey:
                try:
                    import readline
                    readline.set_completer(self.old_completer)
                except ImportError:
                    pass  
            #if not (stop==None):
                #print 'really stop'
                #return stop
                #if stop==_STOP_AND_EXIT:
                #    _STOP_AND_EXIT
                #    print 'I really should stop'
                #    sys.exit()
            return stop

    def do_EOF(self, arg):
        return True
    do_eof = do_EOF
               
    def showParam(self, param):
        any_shown = False
        param = param.strip().lower()
        for p in self.settable:
            if p.startswith(param):
                val = getattr(self, p)
                self.stdout.write('%s: %s\n' % (p, str(getattr(self, p))))
                any_shown = True
        if not any_shown:
            self.perror("Parameter '%s' not supported (type 'show' for list of parameters)." % param)
                
    def do_quit(self, arg):
        return self._STOP_AND_EXIT
    do_exit = do_quit
    do_q = do_quit
    
    def do_show(self, arg):
        '''Shows value of a parameter.'''
        if arg.strip():
            self.showParam(arg)
        else:
            for param in self.settable:
                self.showParam(param)
    
    def do_set(self, arg):
        '''Sets a cmd2 parameter.  Accepts abbreviated parameter names so long as there is no ambiguity.
           Call without arguments for a list of settable parameters with their values.'''
        try:
            paramName, val = arg.split(None, 1)
            paramName = paramName.strip().lower()
            hits = [paramName in p for p in self.settable]
            if hits.count(True) == 1:
                paramName = self.settable[hits.index(True)]
                currentVal = getattr(self, paramName)
                if (val[0] == val[-1]) and val[0] in ("'", '"'):
                    val = val[1:-1]
                else:                
                    val = cast(currentVal, val)
                setattr(self, paramName, val)
                self.stdout.write('%s - was: %s\nnow: %s\n' % (paramName, currentVal, val))
                if currentVal != val:
                    try:
                        onchange_hook = getattr(self, '_onchange_%s' % paramName)
                        onchange_hook(old=currentVal, new=val)
                    except AttributeError:
                        pass
            else:
                self.do_show(paramName)
        except (ValueError, AttributeError, NotSettableError), e:
            self.do_show(arg)
                
    def do_pause(self, arg):
        'Displays the specified text then waits for the user to press RETURN.'
        raw_input(arg + '\n')
        
    def do_shell(self, arg):
        'execute a command as if at the OS prompt.'
        os.system(arg)
                
    def do_py(self, arg):  
        '''
        py <command>: Executes a Python command.
        py: Enters interactive Python mode.
        End with `Ctrl-D` (Unix) / `Ctrl-Z` (Windows), `quit()`, 'exit()`.
        Non-python commands can be issued with `cmd("your command")`.
        '''
        arg = arg.parsed.raw[2:].strip()
        if arg.strip():
            interp = InteractiveInterpreter(locals=self.pystate)
            interp.runcode(arg)
        else:
            interp = MyInteractiveConsole(locals=self.pystate)
            def quit():
                raise EmbeddedConsoleExit
            def onecmd(arg):
                return self.onecmd(arg + '\n')
            self.pystate['quit'] = quit
            self.pystate['exit'] = quit
            self.pystate[self.nonpythoncommand] = onecmd
            try:
                cprt = 'Type "help", "copyright", "credits" or "license" for more information.'        
                keepstate = Statekeeper(sys, ('stdin','stdout'))
                sys.stdout = self.stdout
                sys.stdin = self.stdin
                interp.interact(banner= "Python %s on %s\n%s\n(%s)\n%s" %
                       (sys.version, sys.platform, cprt, self.__class__.__name__, self.do_py.__doc__))
            except EmbeddedConsoleExit:
                pass
            keepstate.restore()
            
    def do_history(self, arg):
        """history [arg]: lists past commands issued
        
        no arg -> list all
        arg is integer -> list one history item, by index
        arg is string -> string search
        arg is /enclosed in forward-slashes/ -> regular expression search
        """
        if arg:
            history = self.history.get(arg)
        else:
            history = self.history
        for hi in history:
            self.stdout.write(hi.pr())
    def last_matching(self, arg):
        try:
            if arg:
                return self.history.get(arg)[-1]
            else:
                return self.history[-1]
        except IndexError:
            return None        
    def do_list(self, arg):
        """list [arg]: lists last command issued
        
        no arg -> list absolute last
        arg is integer -> list one history item, by index
        - arg, arg - (integer) -> list up to or after #arg
        arg is string -> list last command matching string search
        arg is /enclosed in forward-slashes/ -> regular expression search
        """
        try:
            self.stdout.write(self.last_matching(arg).pr())
        except:
            pass
    do_hi = do_history
    do_l = do_list
    do_li = do_list
        
    def do_ed(self, arg):
        """ed: edit most recent command in text editor
        ed [N]: edit numbered command from history
        ed [filename]: edit specified file name
        
        commands are run after editor is closed.
        "set edit (program-name)" or set  EDITOR environment variable
        to control which editing program is used."""
        if not self.editor:
            self.perror("Please use 'set editor' to specify your text editing program of choice.")
            return
        filename = self.default_file_name
        if arg:
            try:
                buffer = self.last_matching(int(arg))
            except ValueError:
                filename = arg
                buffer = ''
        else:
            buffer = self.history[-1]

        if buffer:
            f = open(os.path.expanduser(filename), 'w')
            f.write(buffer or '')
            f.close()        
                
        os.system('%s %s' % (self.editor, filename))
        self.do__load(filename)
    do_edit = do_ed
    
    saveparser = (pyparsing.Optional(pyparsing.Word(pyparsing.nums)^'*')("idx") + 
                  pyparsing.Optional(pyparsing.Word(legalChars + '/\\'))("fname") +
                  pyparsing.stringEnd)    
    def do_save(self, arg):
        """`save [N] [filename.ext]`
        Saves command from history to file.
        N => Number of command (from history), or `*`; 
             most recent command if omitted"""

        try:
            args = self.saveparser.parseString(arg)
        except pyparsing.ParseException:
            self.perror(self.do_save.__doc__)
            return
        fname = args.fname or self.default_file_name
        if args.idx == '*':
            saveme = '\n\n'.join(self.history[:])
        elif args.idx:
            saveme = self.history[int(args.idx)-1]
        else:
            saveme = self.history[-1]
        try:
            f = open(os.path.expanduser(fname), 'w')
            f.write(saveme)
            f.close()
            self.pfeedback('Saved to %s' % (fname))
        except Exception, e:
            self.perror('Error saving %s: %s' % (fname, str(e)))
            
    def read_file_or_url(self, fname):
        if isinstance(fname, file):
            target = open(fname, 'r')
        else:
            match = self.urlre.match(fname)
            if match:
                target = urllib.urlopen(match.group(1))
            else:
                fname = os.path.expanduser(fname)
                try:
                    result = open(os.path.expanduser(fname), 'r')
                except IOError, e:                    
                    result = open('%s.%s' % (os.path.expanduser(fname), 
                                             self.defaultExtension), 'r')
        return result
        
    def do__relative_load(self, arg=None):
        '''
        Runs commands in script at file or URL; if this is called from within an
        already-running script, the filename will be interpreted relative to the 
        already-running script's directory.'''
        if arg:
            arg = arg.split(None, 1)
            targetname, args = arg[0], (arg[1:] or [''])[0]
            targetname = os.path.join(self.current_script_dir or '', targetname)
            self.do__load('%s %s' % (targetname, args))
    
    urlre = re.compile('(https?://[-\\w\\./]+)')
    def do_load(self, arg=None):           
        """Runs script of command(s) from a file or URL."""
        if arg is None:
            targetname = self.default_file_name
        else:
            arg = arg.split(None, 1)
            targetname, args = arg[0], (arg[1:] or [''])[0].strip()
        #try:
        #    target = self.read_file_or_url(targetname)
        #except IOError, e:
        #    self.perror('Problem accessing script from %s: \n%s' % (targetname, e))
        #    return
       
        
        
        keepstate = Statekeeper(self, ('stdin','use_rawinput','prompt',
                                       'continuation_prompt','current_script_dir'))
        flag=True
        while flag:
            try:
                targetstr,i,mark=sequenceloader.driver(targetname)
                target=open('tmp','w')
                target.write(targetstr)
                #target.write(r'\n')
                target.close()
                target=open('tmp','r')
            except IOError, e:
                self.perror('Problem accessing script from %s: \n%s' % (targetname, e))
                return
            self.stdin = target
            #print 'targetstr',targetstr
            self.use_rawinput = False
            self.prompt = self.continuation_prompt = ''
            self.current_script_dir = os.path.split(targetname)[0]
            stop = self.cmdloop()
            self.stdin.close()
            if i==mark:
                break
        #print 'restoring'
        keepstate.restore()
        self.lastcmd = ''
        return (stop == self._STOP_AND_EXIT) and self._STOP_AND_EXIT    
    do__load = do_load  # avoid an unfortunate legacy use of do_load from sqlpython
    
    def do_run(self, arg):
        """run [arg]: re-runs an earlier command
        
        no arg -> run most recent command
        arg is integer -> run one history item, by index
        arg is string -> run most recent command by string search
        arg is /enclosed in forward-slashes/ -> run most recent by regex
        """        
        'run [N]: runs the SQL that was run N commands ago'
        runme = self.last_matching(arg)
        self.pfeedback(runme)
        if runme:
            runme = self.precmd(runme)
            stop = self.onecmd(runme)
            stop = self.postcmd(stop, runme)
    do_r = do_run        
            
    def fileimport(self, statement, source):
        try:
            f = open(os.path.expanduser(source))
        except IOError:
            self.stdout.write("Couldn't read from file %s\n" % source)
            return ''
        data = f.read()
        f.close()
        return data
            
class HistoryItem(str):
    def __init__(self, instr):
        str.__init__(self)
        self.lowercase = self.lower()
        self.idx = None
    def pr(self):
        return '-------------------------[%d]\n%s\n' % (self.idx, str(self))
        
class History(list):
    rangeFrom = re.compile(r'^([\d])+\s*\-$')
    def append(self, new):
        new = HistoryItem(new)
        list.append(self, new)
        new.idx = len(self)
    def extend(self, new):
        for n in new:
            self.append(n)
    def get(self, getme):
        try:
            getme = int(getme)
            if getme < 0:
                return self[:(-1 * getme)]
            else:
                return [self[getme-1]]
        except IndexError:
            return []
        except (ValueError, TypeError):
            getme = getme.strip()
            mtch = self.rangeFrom.search(getme)
            if mtch:
                return self[(int(mtch.group(1))-1):]
            if getme.startswith(r'/') and getme.endswith(r'/'):
                finder = re.compile(getme[1:-1], re.DOTALL | re.MULTILINE | re.IGNORECASE)
                def isin(hi):
                    return finder.search(hi)
            else:
                def isin(hi):
                    return (getme.lower() in hi.lowercase)
            return [itm for itm in self if isin(itm)]

class NotSettableError(Exception):
    pass
        
def cast(current, new):
    """Tries to force a new value into the same type as the current."""
    typ = type(current)
    if typ == bool:
        try:
            return bool(int(new))
        except ValueError, TypeError:
            pass
        try:
            new = new.lower()    
        except:
            pass
        if (new=='on') or (new[0] in ('y','t')):
            return True
        if (new=='off') or (new[0] in ('n','f')):
            return False
    else:
        try:
            return typ(new)
        except:
            pass
    print "Problem setting parameter (now %s) to %s; incorrect type?" % (current, new)
    return current
        
class Statekeeper(object):
    def __init__(self, obj, attribs):
        self.obj = obj
        self.attribs = attribs
        self.save()
    def save(self):
        for attrib in self.attribs:
            setattr(self, attrib, getattr(self.obj, attrib))
    def restore(self):
        for attrib in self.attribs:
            setattr(self.obj, attrib, getattr(self, attrib))        

class Borg(object):
    '''All instances of any Borg subclass will share state.
    from Python Cookbook, 2nd Ed., recipe 6.16'''
    _shared_state = {}
    def __new__(cls, *a, **k):
        obj = object.__new__(cls, *a, **k)
        obj.__dict__ = cls._shared_state
        return obj
    
class OutputTrap(Borg):
    '''Instantiate an OutputTrap to divert/capture ALL stdout output.  For use in unit testing.
    Call `tearDown()` to return to normal output.'''
    def __init__(self):
        self.old_stdout = sys.stdout
        self.trap = tempfile.TemporaryFile()
        sys.stdout = self.trap
    def read(self):
        self.trap.seek(0)
        result = self.trap.read()
        self.trap.truncate(0)
        return result.strip('\x00')        
    def tearDown(self):
        sys.stdout = self.old_stdout

       
class Cmd2TestCase(unittest.TestCase):
    '''Subclass this, setting CmdApp, to make a unittest.TestCase class
       that will execute the commands in a transcript file and expect the results shown.
       See example.py'''
    CmdApp = None
    def fetchTranscripts(self):
        self.transcripts = {}
        for fileset in self.CmdApp.testfiles:
            for fname in glob.glob(fileset):
                tfile = open(fname)
                self.transcripts[fname] = iter(tfile.readlines())
                tfile.close()
        if not len(self.transcripts):
            raise StandardError, "No test files found - nothing to test."
    def setUp(self):
        if self.CmdApp:
            self.outputTrap = OutputTrap()
            self.cmdapp = self.CmdApp()
            self.fetchTranscripts()
    def testall(self):
        if self.CmdApp:
            its = sorted(self.transcripts.items())
            for (fname, transcript) in its:
                self._test_transcript(fname, transcript)
    regexPattern = pyparsing.QuotedString(quoteChar=r'/', escChar='\\', multiline=True, unquoteResults=True)
    regexPattern.ignore(pyparsing.cStyleComment)
    notRegexPattern = pyparsing.Word(pyparsing.printables)
    notRegexPattern.setParseAction(lambda t: re.escape(t[0]))
    expectationParser = regexPattern | notRegexPattern
    endStrippingRegex = re.compile(r'[ \t]*\n')
    def _test_transcript(self, fname, transcript):
        lineNum = 0
        try:
            line = transcript.next()
            while True:
                while not line.startswith(self.cmdapp.prompt):
                    line = transcript.next()
                command = [line[len(self.cmdapp.prompt):]]
                line = transcript.next()
                while line.startswith(self.cmdapp.continuation_prompt):
                    command.append(line[len(self.cmdapp.continuation_prompt):])
                    line = transcript.next()
                command = ''.join(command)
                self.cmdapp.onecmd(command)
                result = self.outputTrap.read().strip()
                if line.startswith(self.cmdapp.prompt):
                    message = '\nFile %s, line %d\nCommand was:\n%s\nExpected: (nothing)\nGot:\n%s\n'%\
                        (fname, lineNum, command, result)     
                    self.assert_(not(result.strip()), message)
                    continue
                expected = []
                while not line.startswith(self.cmdapp.prompt):
                    expected.append(line)
                    line = transcript.next()
                expected = ''.join(expected).strip()
                message = '\nFile %s, line %d\nCommand was:\n%s\nExpected:\n%s\nGot:\n%s\n'%\
                    (fname, lineNum, command, expected, result)      
                expected = self.expectationParser.transformString(expected)
                expected = self.endStrippingRegex.sub('\s*\n', expected)
                self.assert_(re.match(expected, result, re.MULTILINE | re.DOTALL), message)
        except StopIteration:
            pass
    def tearDown(self):
        if self.CmdApp:
            self.outputTrap.tearDown()
        
if __name__ == '__main__':
    doctest.testmod(optionflags = doctest.NORMALIZE_WHITESPACE)
    
'''
To make your application transcript-testable, add text like this to your .py file
(replacing CmdLineApp with your own application class's name).  Then, a cut-and-pasted
version of a successful session with your application, saved as a text file, can serve
as a test for future 

Invoke the test later with `python myapplication.py --test mytranscripttestfile.ext`
Wildcards can be used to test against multiple transcript files.


class TestMyAppCase(Cmd2TestCase):
    CmdApp = CmdLineApp
parser = optparse.OptionParser()
parser.add_option('-t', '--test', dest='test', action="store_true", 
                  help='Test against transcript(s) in FILE (wildcards OK)')
(callopts, callargs) = parser.parse_args()
if callopts.test:
    CmdLineApp.testfiles = callargs
    sys.argv = [sys.argv[0]] # the --test argument upsets unittest.main()
    unittest.main()
else:
    CmdLineApp().cmdloop()
'''