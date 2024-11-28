IMAGE=atlantis_docker
USER=$(id -u)
GROUP=$(id -g)

run: build
	 docker run -v /Users/ful083/AtlantisRepository/AtlantisCurrentTrunk/example:/app/model --user $(id -u):$(id -g) atlantis_docker

#REGISTRY=docker-registry.it.csiro.au
build:
	docker build -t $(IMAGE) .

clean:
	docker system prune -f
