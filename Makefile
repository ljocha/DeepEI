
TAG:="ljocha/deepei"
UID:=$(shell id -u)

build:
	docker build -t ${TAG} .


push:
	docker push ${TAG}


run:
	docker run -u ${UID} -ti -v $$PWD:/work -w /work -e PYTHONPATH=/work ${TAG} bash
