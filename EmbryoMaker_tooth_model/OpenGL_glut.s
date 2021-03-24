	.file	"OpenGL_glut.f90"
	.comm	glut_bitmap_8_by_13,8,8
	.comm	glut_bitmap_9_by_15,8,8
	.comm	glut_bitmap_helvetica_10,8,8
	.comm	glut_bitmap_helvetica_12,8,8
	.comm	glut_bitmap_helvetica_18,8,8
	.comm	glut_bitmap_times_roman_10,8,8
	.comm	glut_bitmap_times_roman_24,8,8
	.comm	glut_stroke_mono_roman,8,8
	.comm	glut_stroke_roman,8,8
	.text
	.globl	glut_null_func_
	.type	glut_null_func_, @function
glut_null_func_:
.LFB0:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	nop
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE0:
	.size	glut_null_func_, .-glut_null_func_
	.ident	"GCC: (Ubuntu 5.4.0-6ubuntu1~16.04.5) 5.4.0 20160609"
	.section	.note.GNU-stack,"",@progbits
