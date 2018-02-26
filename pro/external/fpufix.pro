;+
; NAME:
;       FPUFIX
;
; PURPOSE:
;
;       This is a utility routine to examine a variable and fix problems
;       that will create floating point underflow errors.
;
; AUTHOR:
;
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY:
;
;       Utilities
;
; CALLING SEQUENCE:
;
;       fixedData = FPUFIX(data)
;
; ARGUMENTS:
;
;       data :         A numerical variable to be checked for values that will cause
;                      floating point underflow errors. Suspect values are set to 0.
;
; KEYWORDS:
;
;       None.

; RETURN VALUE:
;
;       fixedData:    The output is the same as the input, except that any values that
;                     will cause subsequent floating point underflow errors are set to 0.
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLES:
;
;       data = FPTFIX(data)
;
; RESTRICTIONS:
;
;     None.
;
; MODIFICATION HISTORY:
;
;       Written by David W. Fanning, from Mati Meron's example FPU_FIX. Mati's
;          program is more robust that this (ftp://cars3.uchicago.edu/midl/),
;          but this serves my needs and doesn't require other programs from
;          Mati's library.  24 February 2006.
;-
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 2006 Fanning Software Consulting
;
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
;
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
;
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;
; 3. This notice may not be removed or altered from any source distribution.
;
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;
;###########################################################################
FUNCTION FPUFIX, data

   ; Return to caller on error after setting !Except.
   On_Error, 2
   Catch, theError
   IF theError NE 0 THEN BEGIN
      Catch, /Cancel
      Message, !Error_State.Msg
      IF N_Elements(except) THEN !EXCEPT = except
      RETURN, data
   ENDIF

   ; Only want to deal with numerical data types.
   ; Return all other kinds.
   dataType = Size(data, /Type)
   nogoodtypes = [0,7,8,10,11]
   void = Where(nogoodTypes EQ dataType, count)
   IF count GT 0 THEN RETURN, data

   ; Floating underflow error we are trying to fix.
   fpu_error = 32

   ; Save current !EXCEPT. Don't report exceptions here.
   except = !EXCEPT
   !EXCEPT = 1

   ; Clear math error status.
   void = Check_Math()

   ; Do something with the data that will cause floating underflow errors.
   void = Min(data, /NAN)

   ; Check the math error status now.
   check = Check_Math()

   ; If this is a floating underflow error, then fix it.
   IF check EQ fpu_error THEN BEGIN
      info = MaChar(DOUBLE=(dataType EQ 5 OR dataType EQ 9))
      indices = Where(Abs(data) LT info.xmin, count)
      IF count GT 0 THEN data[indices] = 0
   ENDIF

   ; Clean up.
   !EXCEPT = except
   void = Check_Math()

   ; Return the repaired data.
   RETURN, data

END
