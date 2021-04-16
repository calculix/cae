def forceFocus(wnd):
	SPI_GETFOREGROUNDLOCKTIMEOUT = 0x2000
	SPI_SETFOREGROUNDLOCKTIMEOUT = 0x2001

	SW_RESTORE = 9
	SPIF_SENDCHANGE = 2
	
	import ctypes
	IsIconic = ctypes.windll.user32.IsIconic
	ShowWindow = ctypes.windll.user32.ShowWindow
	GetForegroundWindow = ctypes.windll.user32.GetForegroundWindow
	GetWindowThreadProcessId = ctypes.windll.user32.GetWindowThreadProcessId
	BringWindowToTop = ctypes.windll.user32.BringWindowToTop
	AttachThreadInput = ctypes.windll.user32.AttachThreadInput
	SetForegroundWindow = ctypes.windll.user32.SetForegroundWindow
	
	if IsIconic(wnd):
		ShowWindow(wnd, SW_RESTORE)

	if GetForegroundWindow() == wnd:
		return True
	
	ForegroundThreadID = GetWindowThreadProcessId(GetForegroundWindow(), None)
	ThisThreadID = GetWindowThreadProcessId(wnd, None)
	if AttachThreadInput(ThisThreadID, ForegroundThreadID, True):
		BringWindowToTop(wnd)
		SetForegroundWindow(wnd)
		AttachThreadInput(ThisThreadID, ForegroundThreadID, False)
		if GetForegroundWindow() == wnd:
			return True

	timeout = ctypes.c_int()
	zero = ctypes.c_int(0)
	# SystemParametersInfo(SPI_GETFOREGROUNDLOCKTIMEOUT, 0, ctypes.byref(timeout), 0)
	# SystemParametersInfo(SPI_SETFOREGROUNDLOCKTIMEOUT, 0, ctypes.byref(zero), SPIF_SENDCHANGE)
	BringWindowToTop(wnd)
	SetForegroundWindow(wnd)
	# SystemParametersInfo(SPI_SETFOREGROUNDLOCKTIMEOUT, 0, ctypes.byref(timeout), SPIF_SENDCHANGE); 
	if GetForegroundWindow() == wnd:
		return True
		
	return False