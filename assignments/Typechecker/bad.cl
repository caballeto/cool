class C {
	a : Int;
	b : Bar;
	init(x : Dar, y : Bool) : C {
           {
		a <- x;
		b <- y;
		self;
           }
	};
};

class Int {
	cmain():C {
	 {
	  (new C).init(1,1);
	  (new C).init(1,true,3);
	  (new C).iinit(1,true);
	  (new C);
	 }
	};
};

