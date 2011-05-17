dojo.provide("plugins.core.ExpandoPane");

dojo.require("dojox.layout.ExpandoPane");

dojo.declare("plugins.core.ExpandoPane",
	[ dojox.layout.ExpandoPane ],
{

	// summary: An adaptation of dojox.layout.ExpandoPane to allow the middle
	// 			pane to be shown/hidden, with corresponding adjustments to the
	//			width of the right pane
	//
	width : null,
	minWidth : 15,
	height : null,
	minHeight : 15,
	expand : false,

	//postCreate : function ()
	//{
	//
	//	this.inherited(arguments);
	//	
	//	//this.expand = arguments.expand;		
	//
	//},

	resize : function(node, newSize){

		if ( newSize )
		{
		}

		// summary: we aren't a layout widget, but need to act like one:
		var size = dojo.marginBox(this.domNode);
		var h = size.h - this._titleHeight;
		dojo.style(this.containerNode, "height", h + "px");

		if ( newSize && newSize.w )
		{
			dojo.style(this.domNode, "width", newSize.w + "px");
		}

		//this.inherited(arguments);

	},


	_showEnd : function()
	{

		// summary: Common animation onEnd code - "unclose"	
		dojo.style(this.cwrapper, { opacity: 0, visibility:"visible" });		
		dojo.fadeIn({ node:this.cwrapper, duration:227 }).play(1);
		dojo.removeClass(this.domNode, "dojoxExpandoClosed");


		if (this.region)
		{
			switch (this.region)
			{
				case "left" : case "center" : case "right" :

					this.domNode.style.width = this.width + "px";
					break;
				case "top" : case "bottom" :
					this.domNode.style.height = this.height + "px";
					break;
			}
		}

		setTimeout(dojo.hitch(this._container, "layout", this.region, this, "show"), 15);
	},


	_hideEnd : function(){

		if (this.region)
		{
			switch (this.region)
			{
				case "left" : case "center" : case "right" :
					dojo.style(this.domNode, "width", this.minWidth + "px");
					//this.domNode.style.width = this.minWidth + "px";
										break;
				case "top" : case "bottom" :
					dojo.style(this.domNode, "width", this.minHeight + "px");
					//this.domNode.style.height = this.minHeight + "px";
					break;
			}
		}

		setTimeout(dojo.hitch(this._container, "layout", this.region, this, "hide" ), 15);
	}

});
