	.file	"OpenGL_gl.f90"
	.text
	.globl	__opengl_gl_MOD_cstring
	.type	__opengl_gl_MOD_cstring, @function
__opengl_gl_MOD_cstring:
.LFB0:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -40(%rbp)
	movl	%esi, -44(%rbp)
	movq	%rdx, -56(%rbp)
	movl	%ecx, -48(%rbp)
	movl	-48(%rbp), %eax
	cltq
	movq	%rax, %r10
	movl	$0, %r11d
	movq	-40(%rbp), %rax
	movq	24(%rax), %rax
	testq	%rax, %rax
	jne	.L3
	movl	$1, %eax
.L3:
	movq	-40(%rbp), %rdx
	movq	(%rdx), %rdx
	movq	%rdx, -16(%rbp)
	movl	-48(%rbp), %edx
	addl	$1, %edx
	movslq	%edx, %rdx
	movq	%rdx, -24(%rbp)
	movq	%rax, %rdx
	imulq	-24(%rbp), %rdx
	movq	%rax, %rcx
	negq	%rcx
	movq	%rdx, %r8
	movl	$0, %r9d
	movl	-48(%rbp), %edx
	movl	$1, -4(%rbp)
	cmpl	%edx, -4(%rbp)
	jg	.L4
.L5:
	movl	-4(%rbp), %esi
	movslq	%esi, %rsi
	imulq	%rax, %rsi
	leaq	(%rsi,%rcx), %r8
	movq	-56(%rbp), %rdi
	movl	-4(%rbp), %esi
	subl	$1, %esi
	movslq	%esi, %rsi
	movzbl	(%rdi,%rsi), %edi
	movq	-16(%rbp), %rsi
	movb	%dil, (%rsi,%r8)
	cmpl	%edx, -4(%rbp)
	sete	%sil
	movzbl	%sil, %esi
	addl	$1, -4(%rbp)
	testl	%esi, %esi
	jne	.L4
	jmp	.L5
.L4:
	movl	-48(%rbp), %edx
	addl	$1, %edx
	movslq	%edx, %rdx
	imulq	%rdx, %rax
	leaq	(%rax,%rcx), %rdx
	movq	-16(%rbp), %rax
	movb	$0, (%rax,%rdx)
	nop
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE0:
	.size	__opengl_gl_MOD_cstring, .-__opengl_gl_MOD_cstring
	.ident	"GCC: (Ubuntu 5.4.0-6ubuntu1~16.04.5) 5.4.0 20160609"
	.section	.note.GNU-stack,"",@progbits
