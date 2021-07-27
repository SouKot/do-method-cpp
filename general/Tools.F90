module m_Tools

implicit none

contains

subroutine Error(string,file,line)

character(len=*) :: string
character(len=*) :: file
integer :: line

write(*,*) 'Error: '//string
write(*,*) '(in file '//file//', line ',line,')'
stop 'error encountered.'

end subroutine Error

subroutine Warning(string,file,line)

character(len=*) :: string
character(len=*) :: file
integer :: line

write(*,*) 'Warning: '//string
write(*,*) '(in file '//file//', line ',line,')'

end subroutine Warning

function toString(int)

implicit none

integer :: int

character(len=4) :: toString
integer :: i

write(toString,'(1I4)') int
do i=1,4
if (toString(i:i)==' ') toString(i:i)='_'
end do
return
end function toString

end module m_Tools
