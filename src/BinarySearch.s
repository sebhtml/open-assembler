.text
.globl BinarySearch
BinarySearch:
	movl	%edx, %r8d
	subl	$1, %r8d
	xorl %ecx,%ecx
	jmp	while
lower:
	leal	-1(%rdx), %r8d
	cmpl	%ecx, %r8d
	jl	notFound
while:
	leal	(%r8,%rcx), %eax
	movl	%eax, %edx
	shrl	$31, %edx
	addl	%eax, %edx
	sarl	%edx
	movslq	%edx,%rax
	movl	%edx, %r9d
	cmpq	%rsi, (%rdi,%rax,8)
	ja	lower
	jae	found
	leal	1(%rdx), %ecx
	cmpl	%ecx, %r8d
	jge	while
notFound:
	movl	$-1, %r9d
found:
	movl	%r9d, %eax
	ret
