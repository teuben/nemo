	interface to subroutine C_qdevice[c](dev)
	integer dev
	end

	interface to subroutine C_unqdevice[c](dev)
	integer dev
	end

	interface to function C_qread[c](val)
	integer *2 val
	integer C_qread
	end

	interface to function C_qtest[c]()
	integer C_qtest
	end

	interface to subroutine C_qreset[c]()
	end

	interface to function C_isqueued[c]()
	integer C_isqueued
	end

	interface to subroutine C_qenter[c](dev, val)
	integer	dev
	integer *2 val
	end

	subroutine qdevice(dev)
	integer dev
	call C_qdevice(dev)
	end

	subroutine qdevic(dev)
	integer dev
	call C_qdevice(dev)
	end

	subroutine unqdevice(dev)
	integer dev
	call C_unqdevice(dev)
	end

	subroutine unqdev(dev)
	integer dev
	call C_unqdevice(dev)
	end

	integer function qtest()
	integer C_qtest
	qtest = C_qtest()
	return
	end

	subroutine qreset()
	call C_qreset
	return
	end

	integer function isqueued()
	integer C_isqueued
	isqueued = C_isqueued()
	return
	end

	integer function isqueu()
	integer C_isqueued
	isqueu = C_isqueued()
	return
	end

	integer function qread(val)
	integer*2 val
	integer C_qread
	qread = C_qread(val)
	return
	end

	subroutine qenter(dev, val)
	integer dev, val
	call C_qenter(dev, val)
	return
	end
