#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:     Nolan Dickson
# Last Edit:  2018-02-07 11:24:11
# Contact:    Nolan.Dickson@Canada.ca

import sys
import time
import itertools
import threading

# TODO is breaking now: stuttering (when too many update calls?)
# TODO also breaking now: model init bar doesnt reach 100%

# TODO stop when there is an error
# TODO a pause and restart method


class Load_message(object):
    '''
    load = pyloader.Load_message('My loading message')
    ...
    load.update('My second loading message')
    ...
    load.stop()
    print('time taken:', load.delta_t')
    '''

    def __init__(self, mssg):
        self._load_iter = [mssg + d for d in ['.', '..', '...', '..']]

        self.done = False
        self.time0 = time.time()

        self._thread = threading.Thread(target=self._write)
        self._thread.daemon = True
        self._thread.start()

    def _write(self):
        for c in itertools.cycle(range(len(self._load_iter))):
            sys.stdout.write(self._load_iter[c] + '    \r')
            sys.stdout.flush()
            time.sleep(0.15)

            if self.done:
                break

    def update(self, mssg):
        sys.stdout.write('\r' + ' ' * len(self._load_iter[0]) + '    \r')
        self._load_iter = [mssg + d for d in ['.', '..', '...', '..']]

    def pause(self, hide=True):
        if hide:
            sys.stdout.write('\r' + ' ' * len(self._load_iter[0]) + '    \r')
        self.done = True

    def restart(self):
        self.done = False

        self._thread = threading.Thread(target=self._write)
        self._thread.daemon = True
        self._thread.start()

    def stop(self):
        time.sleep(0.15)  # necessary?
        sys.stdout.write('\r' + ' ' * len(self._load_iter[0]) + '    \r')
        self.done = True
        self.delta_t = time.time() - self.time0


class Load_animation(object):
    # TODO more animations

    def _animation_choice(self, animation, symbol='o', mssg=''):
        if animation == 'bounce':
            # TODO dynamic for size?
            self._load_iter = ('[{o}------] {m}^^[--{o}----] {m}^^[----{o}--] {m}^^[------{o}] {m}' +
                               '^^[----{o}--] {m}^^[--{o}----] {m}').format(o=symbol, m=mssg).split('^^')

        elif animation == 'spin':
            if symbol == 'line':
                self._load_iter = '| {m}^^/ {m}^^- {m}^^\\ {m}^^| {m}^^/ {m}^^- {m}^^\\ {m}'.format(m=mssg).split('^^')
            elif symbol == 'balloon':
                self._load_iter = '. {m}^^o {m}^^O {m}^^@ {m}^^* {m}'.format(m=mssg).split('^^')

    def __init__(self, animation, symbol='o', mssg=''):

        self._animation_choice(animation, symbol, mssg)

        self.done = False
        self.time0 = time.time()

        self._thread = threading.Thread(target=self._write)
        self._thread.daemon = True
        self._thread.start()

    def _write(self):
        for c in itertools.cycle(range(len(self._load_iter))):
            sys.stdout.write(self._load_iter[c] + '    \r')
            sys.stdout.flush()
            time.sleep(0.15)

            if self.done:
                break

    def update(self, animation='', symbol='', mssg=''):
        # TODO if any left blank?
        sys.stdout.write('\r' + ' ' * len(self._load_iter[0]) + '    \r')
        self._animation_choice(animation, symbol, mssg)

    def stop(self):
        time.sleep(0.15)  # necessary?
        sys.stdout.write('\r' + ' ' * len(self._load_iter[0]) + '    \r')
        self.done = True
        self.delta_t = time.time() - self.time0


class Load_bar(object):
    '''
    can override timeint if you desire a faster updating bar, but will obviously be more resource intensive.

    load = pyloader.Load_bar(len(mylist), mssg='Loading')
    ...
    for ind, item in enumerate(mylist):
        load.update(ind)
        ...
    load.stop(persist=True)
    print('time taken:', load.delta_t')
    '''

    def __init__(self, task_length, symbol='=', mssg='', timeint=0.15):

        # 106 to contain the percentage number
        sys.stdout.write('|' + ' ' * 106 + '| {0}\r'.format(mssg))
        sys.stdout.flush()

        self.updated = 0
        self.done = False

        self.timeint = timeint
        self.mssg = mssg
        self.symbol = symbol
        self.task_length = task_length
        self.time0 = time.time()

        self._thread = threading.Thread(target=self._write)
        self._thread.daemon = True
        self._thread.start()

    def _write(self):
        while True:
            if self.updated:
                self.updated -= 1
                sys.stdout.write('|{0}{1:4}%\r'.format(self.symbol * self.prog, self.prog))
                sys.stdout.flush()

            time.sleep(self.timeint)

            if self.done:
                break

    def update(self, progress):
        self.prog = int(progress / (self.task_length * .01) + 2)  # no clue why I need this 2
        self.updated += 1

    def stop(self, persist=False):
        time.sleep(0.15)  # necessary?
        self.done = True
        self.delta_t = time.time() - self.time0
        if persist:
            sys.stdout.write('\n')
        else:
            sys.stdout.write('\r' + ' ' * 108 + ' ' * len(self.mssg) + '\r')
