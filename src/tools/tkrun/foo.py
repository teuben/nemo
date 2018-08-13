#! /bin/env python

def foo(z=0,a=1,b=2,c=3,d=4,e=5,f='test',g='test',h='test',i='test',j='test',k='test',l='test',m='test',n='test',o='test',p='test',q='test',r='test',s='test',t='test',u='test',v='test',w='test',x='test',y='test'):
   """ This function helps barring the foo.
	#> RADIO a= Me,1,3.45
   a is a Radio button to test the Radio Widget
	#> SCALE b=1 0:10:.1
   b is a Scale to test the Scale Widget
	#> CHECK c= mean,sigma,blue,red
   c is a Check button to test the Check Widget
	#> INFILE z=
   z is an Infile option to test the Infile Widget
   #> INFILE d=
   d is an Outfile option to test the Infile Widget
   #> CHECK e= test1,test2,test3
   e is designed to try and break the code
   This function has 3 parameters, a,b,c"""
   
   print 'z=',z
   print 'a=',a
   print 'b=',b
   print 'c=',c
   print 'd=',d
   print 'e=',e
    
def bar(name):
   print "we are going to run name=",name.__name__
   name(0,1,2,3)
   print "and one more time:"
   name(0,c=1,b=2,a=3)
    
if __name__ == '__main__':
   bar(foo)
