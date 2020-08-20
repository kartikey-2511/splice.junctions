
library(LSD)
library(pheatmap)
library(DT)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

pcode="qzuivpeh"

cg=readRDS("cpm/cancer.genes")
z = readRDS("cpm/genes.name.rds")
ref = readRDS("cpm/ref.rds")
z = toupper(z)
colnames(z) = c("ENSEMBL ID", "Gene Name")
n = z
rownames(z) = c(z[,2])

ui <- fluidPage(

	titlePanel("Gene expression for Breast Cancer patients"),

	sidebarLayout(

		sidebarPanel(width=3, 
			textInput("pwd", h3(strong(span("Enter the passcode to view the plots", style="color:#004C99;font-family: 'times'"))), value=""),
			uiOutput("pwdr"),
			br(),
			h4(strong("Note:"), " This site is published and maintained by ", strong("FunGeL Lab, IIT-D")),
			h4("You can request for the passcode at: "),
			h5(strong(span("bb5170057@iitd.ac.in",style="color:blue"))),
			h5(strong("Kartikey Karnatak, ")),
			h5(strong("FunGeL Lab, Biochemical Sciences, ")),
			h5(strong("Indian Institute of Technology, New Delhi, India")),
			br(),
			h5("Please ", strong("search here"), " for valid", strong(" Gene names"), " and  ", strong("ENSEMBL ids.")),
			DTOutput("gene.list")
		),
			
		mainPanel(width=9,

			br(),
			tabsetPanel(type="tab",
				tabPanel(h4(strong("Expression Plots")),
					br(),
					h4("You can enter either gene name or ENSEMBL id. Please refer from the table given on left for correct id/name."),
					br(),
					fluidRow(
						
						column(6,
							textInput("gene1", h4(strong("Enter the first ENSEMBL-Id/Gene-Name")),
									value = "TSPAN6", width='60%')
						),
						column(6,
							textInput("gene2", h4(strong("Enter the second ENSEMBL-Id/Gene-Name")),
									value = "TNMD", width='60%')
						)
					),
					column(6,
						plotOutput("ex.plot1", height=650, width=650)
					),
					column(6,
						plotOutput("ex.plot2", height=650, width=650)
					)
				),
				tabPanel(h4(strong("Correlation Heatmap")),
					br(),
					h4("For proper and distinct viewing, it is advisable to keep the number of genes < 20"),
					fileInput("genes", h4(strong("Upload the file with gene names"))),
					column(6,
						plotOutput("cor.plot1", height=650, width=650)
					),
					column(6,
						plotOutput("cor.plot2", height=650, width=650)
					)
				),
				tabPanel(h4(strong("GO Enrichment")),
					br(),
					column(4,
						selectInput("go_gene", h4(strong("Select the cancer-related gene of interest")),
									choices=c(cg[,2]), selected=cg[1,2]),
						br(),
						selectInput("go_choice", h4(strong("Select the result to be displayed")),
									choices=c("Dot Plot","Enrichment Map","Ridge Plot","GSEA Plot","gse data table"	), selected="Dot Plot"),
						numericInput("go_num", h5(strong("Give the no. of categories to show on plot")),
									value=15),
						numericInput("go_id", h5(strong("Give the gene id for GSEA plot")),
									value=1),
						br(),
						h4(strong(span("Guidelines:", style="color:blue"))),
						h5("1. The genes available for this section are the ", strong("relevant cancer genes"), " published at ", strong("nature.com/articles/s41436-020-0880-8")),
						h5("2. From the 2nd drop-down menu, select the option ", strong("'gse data table'"), " to view the details of ", strong("enriched terms for the selected gene.")),
						h5("3. If the ", strong("no. of categories"), " input is 15, then the 15 ", strong("top most significant enriched terms"), " are displayed on the plots."),
						h5("4. The ", strong("gene id"), " for the GSEA plot should be an ", strong("integer>0 and no more than the number of rows"), " in the 'gse data table'."),
						h5("5. If you would like to get the enrichment data for a gene not currently available, please drop a mail at the contact id provided.")
					),
					column(8,
						uiOutput("err_txt"),
						DTOutput("gse_list"),
						plotOutput("gse_plot", height=800, width=800)
					)
				)
			)
		)
	)
)

server <- function(input, output) {

	gse_obj = reactive({
		s = input$go_gene
		pos = strtoi(cg[s,3])
		gseSet = floor(pos/10)+1
		pos = (pos+1)%%10
		gse = readRDS(paste("gseList.",gseSet,sep=""))[[pos]]
		return(gse)
	})

	output$err_txt = renderUI({
		if (dim(gse_obj())[1]==0)
			tags$h4(strong(span("No enriched terms for this gene for p-value cutoff of 0.05", style="color:red")))
		})

	output$gse_list = renderDT({
		if (input$pwd==pcode) {
			if (input$go_choice=="gse data table") {
				gse = gse_obj()
				if (dim(gse)[1]!=0) {
					gse = as.data.frame(gse)
					rownames(gse) = NULL
					gse = gse[,c(1:5,8)]
					datatable(gse)
				}
			}
		}
	})

	output$gse_plot = renderPlot({
		if (input$pwd==pcode) {
			gse = gse_obj()
			if (dim(gse)[1]!=0) {
				if (input$go_choice=="Dot Plot")
					dotplot(gse, showCategory=input$go_num, split=".sign") + facet_grid(.~.sign)
				else if (input$go_choice=="Enrichment Map")
					emapplot(gse, showCategory=input$go_num)
				else if (input$go_choice=="Ridge Plot")
					ridgeplot(gse, showCategory=input$go_num) + labs(x = "enrichment distribution")
				else if (input$go_choice=="GSEA Plot")
					gseaplot(gse, by = "all", title = gse$Description[input$go_id], geneSetID=input$go_id)
			}
		}
	})

	input_file <- reactive({
    if (is.null(input$genes)) {
      return("")
    }
    else {
    	readLines(input$genes$datapath)
    }
	})

	plot_x = reactive({
		str = toupper(input$gene1)
		if (substr(str,1,4)=="ENSG")
			return(str)
		else
			return(z[str,1])
	})

	plot_y = reactive({
		str = toupper(input$gene2)
		if (substr(str,1,4)=="ENSG")
			return(str)
		else
			return(z[str,1])
	})

	output$pwdr = renderUI({
		if (input$pwd==pcode)
			tags$h5(strong(span("You are good to go :)",style="color:green")))
		else
			tags$h5(strong(span("Incorrect password :( please try again", style="color:red")))
	})

	output$gene.list = renderDT({
		datatable(n)
	})

	output$ex.plot1 = renderPlot({
		s1 = plot_x()
		s2 = plot_y()
		x1 = read.table("cpm/cancer.csv",sep=",",nrows=1,skip=ref[s1,1]-1)
		x2 = read.table("cpm/cancer.csv",sep=",",nrows=1,skip=ref[s2,1]-1)
		v1 = c(unlist(x1))
		v2 = c(unlist(x2))
		if (input$pwd==pcode)
		heatscatter(v1, v2, cexplot=1, cor=TRUE, method='spearman', log="xy", xlab=toupper(input$gene1), ylab=toupper(input$gene2), main="Cancer tissue")
	})

	output$ex.plot2 = renderPlot({
		s1 = plot_x()
		s2 = plot_y()
		x1 = read.table("cpm/normal.csv",sep=",",nrows=1,skip=ref[s1,1]-1)
		x2 = read.table("cpm/normal.csv",sep=",",nrows=1,skip=ref[s2,1]-1)
		v1 = c(unlist(x1))
		v2 = c(unlist(x2))
		if (input$pwd==pcode)
		heatscatter(v1, v2, cexplot=1, cor=TRUE, method='spearman', log="xy", xlab=toupper(input$gene1), ylab=toupper(input$gene2), main="Normal tissue")
	})

	output$cor.plot1 = renderPlot({
		goi = toupper(input_file())
		if (goi[length(goi)] == "")
			goi = goi[1:length(goi)-1]
		l = length(goi)
		if (l > 0) {
			eid = c(z[goi,1])
			m = matrix(nrow=l, ncol=length(read.table("cpm/cancer.csv",sep=",",nrow=1,skip=1)))
			for (i in c(1:l)) {
				s = eid[i]
				x = read.table("cpm/cancer.csv",sep=",",nrows=1,skip=ref[s,1]-1)
				m[i,] = matrix(unlist(x),nrow=1,ncol=length(x))
			}
			m = t(m)
			colnames(m) = c(goi)
			cm = cor(m, method="spearman")
			if (input$pwd==pcode)
			pheatmap(cm, display_numbers = T, main="Cancer tissue")
		}
	})

	output$cor.plot2 = renderPlot({
		goi = toupper(input_file())
		if (goi[length(goi)] == "")
			goi = goi[1:length(goi)-1]
		l = length(goi)
		if (l > 0) {
			eid = matrix(z[goi,1], nrow=l, ncol=1)
			m = matrix(nrow=l, ncol=length(read.table("cpm/normal.csv",sep=",",nrow=1,skip=1)))
			for (i in c(1:l)) {
				s = eid[i]
				x = read.table("cpm/normal.csv",sep=",",nrows=1,skip=ref[s,1]-1)
				m[i,] = matrix(unlist(x),nrow=1,ncol=length(x))
			}
			m = t(m)
			colnames(m) = c(goi)
			cm = cor(m, method="spearman")
			if (input$pwd==pcode)
			pheatmap(cm, display_numbers = T, main="Normal tissue")
		}
	})

}

shinyApp(ui = ui, server = server)
