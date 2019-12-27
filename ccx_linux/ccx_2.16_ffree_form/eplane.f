      !
      !     CalculiX - A 3-dimensional finite element program
      !              Copyright (C) 1998-2019 Guido Dhondt
      !
      !     This program is free software; you can redistribute it and/or
      !     modify it under the terms of the GNU General Public License as
      !     published by the Free Software Foundation(version 2);
      !
      !
      !     This program is distributed in the hope that it will be useful,
      !     but WITHOUT ANY WARRANTY; without even the implied warranty of
      !     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      !     GNU General Public License for more details.
      !
      !     You should have received a copy of the GNU General Public License
      !     along with this program; if not, write to the Free Software
      !     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
      !
      real*8 function eplane(x,n,t)
      !
      implicit none
      !
      real*8 x(*),n(*),t
      !
      eplane= x(1)*n(1)+x(2)*n(2)+x(3)*n(3)+t
      !
      return
      end
      
