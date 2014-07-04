JAVAC=javac
JAVADOC=javadoc
JAR=jar


MAIN_CLASS=VCF2diploid.java
MANIFEST_FILE=DiploidManifest
JAR_FILE=vcf2diploid.jar
CLASS_DIR=classes
DOC_DIR=doc
JAVA_API_URL=http://java.sun.com/j2se/1.5.0/docs/api/

all: prepare compile jar_all

jar_all:
	@rm -f $(JAR_FILE)
	$(JAR) cmf $(MANIFEST_FILE) $(JAR_FILE) -C $(CLASS_DIR) .

compile:
	@echo "Compiling vcf2diploid"
	$(JAVAC) -d $(CLASS_DIR) $(MAIN_CLASS)

prepare:
	@mkdir -p $(CLASS_DIR)

clean:
	@rm -rf $(CLASS_DIR) *.class

