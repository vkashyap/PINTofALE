;+
; NAME:
;       CONVERT_TO_TYPE
;
; PURPOSE:
;
;       Converts its input argument to a specified data type.
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
;       result = Convert_To_Type(input, type)
;
; INPUT_PARAMETERS:
;
;       input:          The input data to be converted.
;       type:           The data type. Accepts values as given by Size(var, /TNAME) or Size(var, /TYPE).
;
; OUTPUT_PARAMETERS:
;
;      result:          The input data is converted to specified data type.
;
; KEYWORDS:
;
;     None.
;
; RESTRICTIONS:
;
;     Data types STRUCT, POINTER, and OBJREF are not allowed.
;
; MODIFICATION HISTORY:
;
;     Written by David W. Fanning, 19 February 2006.
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
FUNCTION Convert_To_Type, input, type

   ; Return to caller on error.
   On_Error, 2

   ; Two positional parameters are required.
   IF N_Params() NE 2 THEN Message, 'Two input parameters (INPUT and TYPE) are required.'

   ; If type is a string, turn it into a number.
   IF Size(type, /TNAME) EQ 'STRING' THEN BEGIN

      type = StrUpCase(type[0])
      CASE type OF
         'BYTE': type = 1
         'INT': type = 2
         'LONG': type = 3
         'FLOAT': type = 4
         'DOUBLE': type = 5
         'COMPLEX': type = 6
         'STRING': type = 7
         'DCOMPLEX': type = 9
         'UNIT': type = 12
         'ULONG': type = 13
         'LONG64': type = 14
         'ULONG64': type = 15
         ELSE: Message, 'Unable to convert input to type: ' + StrUpCase(type)
      ENDCASE

   ENDIF ELSE BEGIN

      ; Only certain kinds of data conversions can occur.
      type = type[0]
      CASE 1 OF
         (type LT 1): Message, 'Unable to convert input to UNDEFINED data type.'
         (type EQ 8): Message, 'Unable to convert input to STRUCTURE data type.'
         (type EQ 10): Message, 'Unable to convert input to POINTER data type.'
         (type EQ 11): Message, 'Unable to convert input to OBJECT data type.'
         (type GT 15): Message, 'Unable to convert undefined data type: ', StrTrim(type) + '.'
         ELSE:
      ENDCASE
   ENDELSE

   ; Do the conversion.
   CASE type OF
      1: output = BYTE(input)
      2: output = FIX(input)
      3: output = LONG(input)
      4: output = FLOAT(input)
      5: output = DOUBLE(input)
      6: output = COMPLEX(input)
      7: output = STRING(input)
      9: output = DCOMPLEX(input)
      12: output = UINT(input)
      13: output = ULONG(input)
      14: output = LONG64(input)
      15: output = ULONG64(input)
   ENDCASE

   RETURN, output

END
