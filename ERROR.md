This error is not written to the log.



WARNING: No slave window.
WARNING: No slave window.
DEBUG: thread_1_read_cgx_stdout_1677673078 STOPPED.
DEBUG: Stopping threads:
thread_1_read_cgx_stdout_1677673078

Traceback (most recent call last):
  File "\\PH-US2\general\PROGRAMMING\WORK\Mirzov\cae\src\actions.py", line 52, in <lambda>
    wf.mw.action_cgx_inp.triggered.connect(lambda: cgx.open_inp(j.inp, len(m.Mesh.nodes)))
  File "S:\CAE\general\PROGRAMMING\WORK\Mirzov\cae\src\gui\cgx.py", line 117, in open_inp
    wf.run_slave(cmd)
  File "S:\CAE\general\PROGRAMMING\WORK\Mirzov\cae\src\gui\window.py", line 148, in run_slave
    self.kill_slave()
  File "S:\CAE\general\PROGRAMMING\WORK\Mirzov\cae\src\gui\window.py", line 187, in kill_slave
    ctypes.windll.user32.SetForegroundWindow(self.sw.info.wid)
AttributeError: 'NoneType' object has no attribute 'wid'
Something went wrong.